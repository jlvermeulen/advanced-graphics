#define _USE_MATH_DEFINES
#define RUSSIAN_ROULETTE_PROBABILITY 0.90f

#define MIN_PATH_LENGTH 2
#define TREE 1 // 0 octtree 1 bvh
//#define MAXLIGHTS 32  

#include "Scene.h"

#include "Intersections.h"
#include "ObjReader.h"
#include <algorithm>
#include <math.h>
#include <random>

Scene::Scene() :
  gen(std::mt19937(std::random_device()())),
  dist(std::uniform_real_distribution<float>(0, 1))
{ }

Scene::~Scene() { }


void Scene::PreRender()
{
	// Instantiate spatial indices
	for (Object* obj : objects)
#if TREE == 1
		obj->ConstructBVHTree();
#elif TREE == 0
		obj->ConstructOctree();
#endif
	// prepare sample light stuff
	unsigned int c = 0;
	for (Object* o : lights)
	{
		for (Triangle* t : o->triangles)
		{
			posv0X[c] = t->Vertices[0].Position.X;
			posv1X[c] = t->Vertices[1].Position.X;
			posv2X[c] = t->Vertices[2].Position.X;
			posv0Y[c] = t->Vertices[0].Position.Y;
			posv1Y[c] = t->Vertices[1].Position.Y;
			posv2Y[c] = t->Vertices[2].Position.Y;
			posv0Z[c] = t->Vertices[0].Position.Z;
			posv1Z[c] = t->Vertices[1].Position.Z;
			posv2Z[c] = t->Vertices[2].Position.Z;

			norv0X[c] = t->Vertices[0].Normal.X;
			norv1X[c] = t->Vertices[1].Normal.X;
			norv2X[c] = t->Vertices[2].Normal.X;
			norv0Y[c] = t->Vertices[0].Normal.Y;
			norv1Y[c] = t->Vertices[1].Normal.Y;
			norv2Y[c] = t->Vertices[2].Normal.Y;
			norv0Z[c] = t->Vertices[0].Normal.Z;
			norv1Z[c] = t->Vertices[1].Normal.Z;
			norv2Z[c] = t->Vertices[2].Normal.Z;

			prev0X[c] = t->v0.X;
			prev0Y[c] = t->v0.Y;
			prev0Z[c] = t->v0.Z;
			prev1X[c] = t->v1.X;
			prev1Y[c] = t->v1.Y;
			prev1Z[c] = t->v1.Z;

			area[c] = t->Area;
			emissionR[c] = o->material.emission.R;
			invDenom[c] = t->invDenom;
			d00[c] = t->d00;
			d01[c] = t->d01;
			d11[c] = t->d11;

			lightTriangles[c] = t;

			c++;
		}
	}

	this->nLight = c;
	this->nLight8 = c / NROFLANES;
	if (c % NROFLANES != 0)
		this->nLight8++;

	// pad with zeroes
	while (c < MAXLIGHTS)
	{
		posv0X[c] = 0.0f;
		posv1X[c] = 0.0f;
		posv2X[c] = 0.0f;
		posv0Y[c] = 0.0f;
		posv1Y[c] = 0.0f;
		posv2Y[c] = 0.0f;
		posv0Z[c] = 0.0f;
		posv1Z[c] = 0.0f;
		posv2Z[c] = 0.0f;

		norv0X[c] = 0.0f;
		norv1X[c] = 0.0f;
		norv2X[c] = 0.0f;
		norv0Y[c] = 0.0f;
		norv1Y[c] = 0.0f;
		norv2Y[c] = 0.0f;
		norv0Z[c] = 0.0f;
		norv1Z[c] = 0.0f;
		norv2Z[c] = 0.0f;

		prev0X[c] = 0.0f;
		prev0Y[c] = 0.0f;
		prev0Z[c] = 0.0f;
		prev1X[c] = 0.0f;
		prev1Y[c] = 0.0f;
		prev1Z[c] = 0.0f;

		area[c] = 0.0f;
		emissionR[c] = 0.0f;
		invDenom[c] = 0.0f;
		d00[c] = 0.0f;
		d01[c] = 0.0f;
		d11[c] = 0.0f;

		lightTriangles[c] = nullptr;

		c++;
	}
}

void Scene::Render(uchar* imageData, int samplesPerPixel, float sigma, bool useDoF)
{
	std::pair<Color3F, float>* samples = new std::pair<Color3F, float>[camera.Width * camera.Height];
	TracePixels(samples, samplesPerPixel, sigma, useDoF);

	#pragma omp parallel for
	for (int x = 0; x < camera.Width; ++x)
		for (int y = 0; y < camera.Height; ++y)
		{
			std::pair<Color3F, float> sample = samples[y * camera.Width + x];
			Color3F color = sample.first / sample.second;
			color.Clip();

			int offset = (y * camera.Width + x) * 4;

			// For each color channel in reversed order (i.e. blue-green-red-alpha)
			imageData[offset] = (uchar) (color.B * 255.0f);
			imageData[offset + 1] = (uchar) (color.G * 255.0f);
			imageData[offset + 2] = (uchar) (color.R * 255.0f);
			imageData[offset + 3] = 255;
		}

	delete [] samples;

	return;
}

void Scene::PostRender()
{

}

//--------------------------------------------------------------------------------
void Scene::TracePixels(std::pair<Color3F, float>* pixelData, int samplesPerPixel, float sigma, bool useDoF)
{
	float tanHalfFovY = tanf(camera.FovY() / 360 * M_PI);
	float tanHalfFovX = tanHalfFovY * camera.Width / camera.Height;

	float left = -tanHalfFovX;
	float right = tanHalfFovX;
	float top = tanHalfFovY;
	float bottom = -tanHalfFovY;

	for (int y = 0; y < camera.Height; ++y)
		for (int x = 0; x < camera.Width; ++x)
			pixelData[y * camera.Width + x] = std::pair<Color3F, float>(Color3F(), 0.0f);

	// Calculate pixel rays
	for (int o = 0; o < 3; ++o) // prevent concurrent access of same array indices
		#pragma omp parallel for schedule(dynamic)
		for (int x = o; x < camera.Width; x += 3)
			for (int y = 0; y < camera.Height; ++y)
				for (int r = 0; r < samplesPerPixel; ++r)
				{
					float rX = dist(gen);
					float rY = dist(gen);

					float a = left + (right - left) * (x + rX) / camera.Width;
					float b = top + (bottom - top) * (y + rY) / camera.Height;

					Vector3F origin = camera.Eye();
					Vector3F direction = Vector3F::Normalise(camera.Focus() + a * camera.Right() + b * camera.Up());

					if (useDoF)
					{
						Vector3F focalPoint = origin + direction * camera.FocalDepth;

						// Get uniformly distributed square [-1,1] x [-1,1]
						float angle = dist(gen) * 2.0f * M_PI;
						float radius = dist(gen);

						float aX = cosf(angle) * radius * camera.Aperture;
						float aY = sinf(angle) * radius * camera.Aperture;

						Vector3F apertureOffset(aX, aY, 0.0f);

						origin += apertureOffset;

						direction *= camera.FocalDepth;
						direction -= apertureOffset;
						direction.Normalise();
					}

					Ray cameraRay(origin, direction);
					Color3F value = TraceRay(cameraRay, true);

					// distribute over neighbouring pixels
					for (int i = -1; i < 2; ++i)
					{
						int xx = x + i;
						if (xx < 0 || xx >= camera.Width)
							continue;

						for (int j = -1; j < 2; ++j)
						{
							int yy = y + j;
							if (yy < 0 || yy >= camera.Height)
								continue;

							float weight = GaussianWeight(i - rX + 0.5f, j - rY + 0.5f, sigma);

							int index = yy * camera.Width + xx;
							std::pair<Color3F, float>& data = pixelData[index];

							data.first += value * weight;
							data.second += weight;
						}
					}
				}
}

//--------------------------------------------------------------------------------
float Scene::GaussianWeight(float dx, float dy, float sigma) const
{
	return exp(-(dx * dx + dy * dy) / (2 * sigma * sigma));
}

//--------------------------------------------------------------------------------
Color3F Scene::TraceRay(const Ray& ray, bool nee)
{
	Material hitMaterial;
	float hitTime;
	Triangle* hitTriangle = FirstHitInfo(ray, hitTime, hitMaterial);
	if (hitTriangle == nullptr)
		return Color3F();

	return ComputeRadiance(ray.Origin + hitTime * ray.Direction, ray.Direction, hitTriangle, hitMaterial, 0, nee);
}

//--------------------------------------------------------------------------------
Color3F Scene::ComputeRadiance(const Vector3F& point, const Vector3F& in, Triangle* triangle, const Material& material, unsigned int depth, bool nee)
{
	const Vector3F& normal = triangle->surfaceNormal(point);
	Color3F c = material.emission;
	if (c.IsSignificant())
		return c;

	if (nee)
		c += DirectIllumination(point, in, normal, material);
	c += IndirectIllumination(point, in, normal, material, depth, nee);

	float mag = c.Magnitude(); // prevent giant weights giving fireflies
	if (mag > 5)
	{
		c /= mag;
		c *= 5;
	}

	return c;
}

Color3F Scene::DirectIllumination(const Vector3F& point, const Vector3F& in, const Vector3F& normal, const Material& material)
{
	if (material.reflType == ReflectionType::specular || material.reflType == ReflectionType::refractive || (material.reflType == ReflectionType::glossy && material.specularExponent > 5))
		return Color3F();

	float lightWeight;
	std::pair<Ray, Triangle*>* sample = SampleLight(point, lightWeight);
	if (sample == nullptr)
		return Color3F();

	Material hitMaterial;
	float hitTime;

	Ray ray = sample->first;
	Triangle* hitTriangle = FirstHitInfo(ray, hitTime, hitMaterial);
	if (hitTriangle == nullptr || hitTriangle != sample->second)
	{
		delete sample;
		return Color3F();
	}

	delete sample;

	Vector3F hitPoint = ray.Origin + hitTime * ray.Direction;
	Vector3F triNormal = hitTriangle->surfaceNormal(hitPoint);
	float weight = lightWeight * std::max(0.0f, Vector3F::Dot(ray.Direction, normal));
	if (material.reflType == ReflectionType::diffuse)
	{
		weight *= 1.0f / M_PI;
	}
	else if (material.reflType == ReflectionType::glossy)
	{
		Vector3F idealReflection = in + 2 * -Vector3F::Dot(normal, in) * normal;
		weight *= pow(std::max(0.0f, Vector3F::Dot(idealReflection, ray.Direction)), 1.0f + material.specularExponent) / pow(M_PI, 1.0f / (1.0f + material.specularExponent));
	}

	return weight * hitMaterial.emission * material.color;
}

Color3F Scene::IndirectIllumination(Vector3F point, const Vector3F& in, const Vector3F& normal, const Material& material, unsigned int depth, bool nee)
{
	if (depth > MIN_PATH_LENGTH && dist(gen) > RUSSIAN_ROULETTE_PROBABILITY)
		return Color3F();

	Ray ray(point, in);
	Color3F value = material.color;

	if (material.reflType == ReflectionType::specular)
	{
		if (Vector3F::Dot(ray.Direction, normal) < 0)
			ray.Reflect(point, normal);
	}
	else if (material.reflType == ReflectionType::glossy)
	{
		// Sample unit sphere
		float phi = 2.0f * M_PI * dist(gen);

		float cosTheta = pow(dist(gen), 1.0f / (1.0f + material.specularExponent)); // importance sampling on phong lobe
		float sinTheta = sqrtf(1.0f - cosTheta * cosTheta);

		float x = sinTheta * cosf(phi);
		float y = sinTheta * sinf(phi);
		float z = cosTheta;

		Vector3F hemi(x, y, z);

		ray.Reflect(point, normal);

		Vector3F w = ray.Direction;

		Vector3F u = Vector3F::Cross(Vector3F(0.00424f, 1, 0.00764f), w);
		u.Normalise();

		Vector3F v = Vector3F::Cross(u, w);

		Vector3F reflDir = hemi.X * u + hemi.Y * v + hemi.Z * w;

		if (Vector3F::Dot(reflDir, normal) < 0) // flip if "behind" normal
			reflDir *= -1;

		ray = Ray(point, reflDir);
		value *= std::max(0.0f, Vector3F::Dot(reflDir, normal));
	}
	else if (material.reflType == ReflectionType::diffuse)
	{
		// cosine importance sampling
		float u1 = dist(gen), u2 = dist(gen);
		float theta = acosf(sqrtf(1.0f - u1));
		float phi = 2.0f * M_PI * u2;

		float xs = sinf(theta) * cosf(phi);
		float ys = cosf(theta);
		float zs = sinf(theta) * sinf(phi);

		Vector3F y(normal.X, normal.Y, normal.Z);
		Vector3F h = y;
		if (fabs(h.X) <= fabs(h.Y) && fabs(h.X) <= fabs(h.Z))
			h.X = 1.0f;
		else if (fabs(h.Y) <= fabs(h.X) && fabs(h.Y) <= fabs(h.Z))
			h.Y = 1.0f;
		else
			h.Z = 1.0f;

		Vector3F x = Vector3F::Normalise(Vector3F::Cross(h, y));
		Vector3F z = Vector3F::Normalise(Vector3F::Cross(x, y));

		Vector3F hemi = xs * x + ys * y + zs * z;
		hemi.Normalise();

		if (Vector3F::Dot(hemi, normal) < 0) // flip if "behind" normal
			hemi *= -1;

		ray = Ray(point, hemi);
	}
	else if (material.reflType == ReflectionType::refractive)
		ray.Refract(point, normal, 1.0f, material.refrIndex);

	Material hitMaterial;
	float hitTime;
	Triangle* hitTriangle = FirstHitInfo(ray, hitTime, hitMaterial);
	if (hitTriangle == nullptr) // didn't hit anything
		return Color3F();

	if (nee && hitMaterial.emission.IsSignificant() &&															// discount hitting a light when using next event estimation
		!(material.reflType == ReflectionType::specular || material.reflType == ReflectionType::refractive ||	// unless we hit a mirror-like surface
			(material.reflType == ReflectionType::glossy && material.specularExponent > 5)))
		return Color3F();

	Vector3F hitPoint = ray.Origin + hitTime * ray.Direction;
	return value * ComputeRadiance(hitPoint, ray.Direction, hitTriangle, hitMaterial, depth + 1, nee) / (depth > MIN_PATH_LENGTH ? RUSSIAN_ROULETTE_PROBABILITY : 1.0f);
}

//--------------------------------------------------------------------------------
Triangle* Scene::FirstHitInfo(const Ray& ray, float& time, Material& mat) const
{
	time = std::numeric_limits<float>::max();
	Triangle* triangle = nullptr;

	for (Object* obj : objects)
	{
		float t = std::numeric_limits<float>::max();
#if TREE == 1
		Triangle* tri = obj->bvhTree->Query(ray, t);
#elif TREE == 0
		Triangle* tri = obj->octree->Query(ray, t);
#endif
		if (tri != nullptr && t < time)
		{
			mat = obj->material;
			time = t;
			triangle = tri;
		}
	}

	return triangle;
}

//--------------------------------------------------------------------------------
std::pair<Ray, Triangle*>* Scene::SampleLight(const Vector3F& hitPoint, float& weight)
{
	weight = 0;

	// create random point within the lightsource
	float u = dist(gen), v = dist(gen);
	if (u + v > 1)
	{
		u = 1 - u;
		v = 1 - v;
	}

	//float flux[MAXLIGHTS];// = new std::pair<float, unsigned int>[lightCount];
	float totalflux = 0;

	// prepare for AVX loop
	const __m256 u8 = _mm256_set1_ps(u);
	const __m256 v8 = _mm256_set1_ps(v);

	const __m256 hitPointX = _mm256_set1_ps(hitPoint.X);
	const __m256 hitPointY = _mm256_set1_ps(hitPoint.Y);
	const __m256 hitPointZ = _mm256_set1_ps(hitPoint.Z);

	__m256 weights = _mm256_setzero_ps();
	__m256 totalflux8 = _mm256_setzero_ps();
	union{ float flux[MAXLIGHTS]; __m256 flux8[MAXLIGHTS / NROFLANES]; };
	union{ float rayX[MAXLIGHTS]; __m256 rayX8[MAXLIGHTS / NROFLANES]; };
	union{ float rayY[MAXLIGHTS]; __m256 rayY8[MAXLIGHTS / NROFLANES]; };
	union{ float rayZ[MAXLIGHTS]; __m256 rayZ8[MAXLIGHTS / NROFLANES]; };

	// loop over all lightsources
	for (unsigned int c = 0; c < nLight8; ++c)
	{
		// Vector3F triPoint = v1 * u + v2 * v + t->Vertices[0].Position;
		const __m256 triPointX = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(prev0X8[c], u8), _mm256_mul_ps(prev1X8[c], v8)), posv0X8[c]);
		const __m256 triPointY = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(prev0Y8[c], u8), _mm256_mul_ps(prev1Y8[c], v8)), posv0Y8[c]);
		const __m256 triPointZ = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(prev0Z8[c], u8), _mm256_mul_ps(prev1Z8[c], v8)), posv0Z8[c]);

		// Vector3F outgoingRay = hitPoint - triPoint;
		__m256 outgoingRayX = _mm256_sub_ps(hitPointX, triPointX);
		__m256 outgoingRayY = _mm256_sub_ps(hitPointY, triPointY);
		__m256 outgoingRayZ = _mm256_sub_ps(hitPointZ, triPointZ);

		// float distance = outgoingRay.LengthSquared();
		const __m256 distance = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(outgoingRayX, outgoingRayX), _mm256_mul_ps(outgoingRayY, outgoingRayY)), _mm256_mul_ps(outgoingRayZ, outgoingRayZ));

		// interpolate surface normal
		// Vector3F v3 = triPoint - t->Vertices[0].Position;
		const __m256 locv3X = _mm256_sub_ps(triPointX,  posv0X8[c]);
		const __m256 locv3Y = _mm256_sub_ps(triPointY,  posv0Y8[c]);
		const __m256 locv3Z = _mm256_sub_ps(triPointZ,  posv0Z8[c]);

		// float d20 = Vector3F::Dot(v3, t->v0);
		const __m256 d20 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(locv3X, prev0X8[c]), _mm256_mul_ps(locv3Y, prev0Y8[c])), _mm256_mul_ps(locv3Z, prev0Z8[c]));

		// float d21 = Vector3F::Dot(v3, t->v1);
		const __m256 d21 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(locv3X, prev1X8[c]), _mm256_mul_ps(locv3Y, prev1Y8[c])), _mm256_mul_ps(locv3Z, prev1Z8[c]));

		// float a2 = (d11 * d20 - d01 * d21) * invDenom;
		const __m256 facY = _mm256_mul_ps(_mm256_sub_ps(_mm256_mul_ps(d118[c], d20), _mm256_mul_ps(d018[c], d21)), invDenom8[c]);

		// float a3 = (d00 * d21 - d01 * d20) * invDenom;
		const __m256 facZ = _mm256_mul_ps(_mm256_sub_ps(_mm256_mul_ps(d008[c], d21), _mm256_mul_ps(d018[c], d20)), invDenom8[c]);

		// float a1 = 1.0f - a2 - a3;
		const __m256 facX = _mm256_sub_ps(_mm256_sub_ps(_mm256_set1_ps(1.0f), facY), facZ);

		// Vector3F normal = factors.X * t->Vertices[0].Normal + factors.Y * t->Vertices[1].Normal + factors.Z * t->Vertices[2].Normal;
		const __m256 normalX = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(facX, norv0X8[c]), _mm256_mul_ps(facY, norv1X8[c])), _mm256_mul_ps(facZ, norv2X8[c]));
		const __m256 normalY = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(facY, norv0Y8[c]), _mm256_mul_ps(facY, norv1Y8[c])), _mm256_mul_ps(facZ, norv2Y8[c]));
		const __m256 normalZ = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(facZ, norv0Z8[c]), _mm256_mul_ps(facY, norv1Z8[c])), _mm256_mul_ps(facZ, norv2Z8[c]));

		// outgoingRay.Normalise();
		const __m256 rdist = _mm256_rsqrt_ps(distance);
		outgoingRayX = _mm256_mul_ps(outgoingRayX, rdist);
		outgoingRayY = _mm256_mul_ps(outgoingRayY, rdist);
		outgoingRayZ = _mm256_mul_ps(outgoingRayZ, rdist);

		//float f = t->Area * std::max(0.0f, Vector3F::Dot(outgoingRay, normal)) / distance;
		const __m256 fdot = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(outgoingRayX, normalX), _mm256_mul_ps(outgoingRayY, normalY)), _mm256_mul_ps(outgoingRayZ, normalZ));
		__m256 f = _mm256_div_ps(_mm256_mul_ps(area8[c], _mm256_max_ps(_mm256_setzero_ps(), fdot)), distance);

		weights = _mm256_add_ps(weights, f); // projected area of all lightsources on hemisphere
		f = _mm256_mul_ps(f, emissionR8[c]); // we assume non-coloured emissions

		flux8[c] = f;
		totalflux8 = _mm256_add_ps(totalflux8, f);

		rayX8[c] = outgoingRayX;
		rayY8[c] = outgoingRayY;
		rayZ8[c] = outgoingRayZ;
	}

	// go from AVX to floats (weights, totalflux)
	// sum parts of totalflux
	for (int i = 0; i < NROFLANES; i++)
	{
		totalflux += totalflux8.m256_f32[i];
		weight += weights.m256_f32[i];
	}

	if (totalflux < 0.001f)
		return nullptr;

	// get the light using the probabilities and a random float between 0 and 1
	float r = dist(gen), invFlux = 1.0f / totalflux;
	int p = -1; // index of chosen triangle
	for (unsigned int i = 0; i < nLight; ++i)
	{
		flux[i] *= invFlux; // dividing all flux values by the total flux gives the probabilities
		if (flux[i] + 0.001f >= r) // add small value to prevent floating point errors
		{
			p = i;
			break;
		}
		else
			r -= flux[i];
	}

	Vector3F v1 = lightTriangles[p]->Vertices[1].Position - lightTriangles[p]->Vertices[0].Position;
	Vector3F v2 = lightTriangles[p]->Vertices[2].Position - lightTriangles[p]->Vertices[0].Position;

	// return ray towards chosen point
	return new std::pair<Ray, Triangle*>(Ray(hitPoint, Vector3F(-rayX[p], -rayY[p], -rayZ[p])), lightTriangles[p]);
}

//--------------------------------------------------------------------------------
void Scene::LoadDefaultScene()
{
	ObjReader reader;
	Clear();
	unsigned int nTriangles, nr = 0;

	// Parse spheres
	reader.parseFile("../models/sphere.obj");
	reader.parseFile("../models/sphere.obj");
	reader.parseFile("../models/sphere.obj");
	reader.parseFile("../models/sphere.obj");
	reader.parseFile("../models/sphere.obj");

	// Parse right
	reader.parseFile("../models/cube.obj");

	// Parse left
	reader.parseFile("../models/cube.obj");

	// Parse back
	reader.parseFile("../models/cube.obj");

	// Parse top
	reader.parseFile("../models/cube.obj");

	// Parse bottom
	reader.parseFile("../models/cube.obj");

	// Parse light
	objects = reader.parseFile("../models/light4.obj");

	// Left sphere
	objects[nr]->material = Material(ReflectionType::glossy, Color3F(0.75f, 0.25f, 0.25f), Color3F(), 1.0f, 100.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position /= 2;

			objects[nr]->triangles[i]->Vertices[j].Position.X += 0.8f;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Left sphere
	objects[nr]->material = Material(ReflectionType::specular, Color3F(1.0f, 1.0f, 1.0f), Color3F(), 1.0f, 0.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position /= 3.5f;

			objects[nr]->triangles[i]->Vertices[j].Position.X += 0.3f;
			objects[nr]->triangles[i]->Vertices[j].Position.Y += 0.8f;
			objects[nr]->triangles[i]->Vertices[j].Position.Z -= 1.0f;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Left sphere
	objects[nr]->material = Material(ReflectionType::refractive, Color3F(1.0f, 1.0f, 1.0f), Color3F(), 1.5f, 0.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position /= 1.5f;

			objects[nr]->triangles[i]->Vertices[j].Position.X -= 0.7f;
			objects[nr]->triangles[i]->Vertices[j].Position.Y -= 5.0f / 6.0f;
			objects[nr]->triangles[i]->Vertices[j].Position.Z -= 0.7f;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Left sphere
	objects[nr]->material = Material(ReflectionType::glossy, Color3F(0.25f, 0.25f, 0.75f), Color3F(), 1.0f, 10.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position /= 3;

			objects[nr]->triangles[i]->Vertices[j].Position.X -= 0.75f;
			objects[nr]->triangles[i]->Vertices[j].Position.Y += 0.75f;
			objects[nr]->triangles[i]->Vertices[j].Position.Z += 2.0f;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Left sphere
	objects[nr]->material = Material(ReflectionType::glossy, Color3F(0.5f, 0.0f, 0.0f), Color3F(), 1.0f, 0.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position /= 4;

			objects[nr]->triangles[i]->Vertices[j].Position.X += 0.5f;
			objects[nr]->triangles[i]->Vertices[j].Position.Y -= 1.25f;
			objects[nr]->triangles[i]->Vertices[j].Position.Z += 0.5f;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Right
	objects[nr]->material = Material(ReflectionType::diffuse, Color3F(1.0f, 0.6f, 0.6f), Color3F(), 1.0f, 0.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position *= 5;
			objects[nr]->triangles[i]->Vertices[j].Position.X += 4;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Left
	objects[nr]->material = Material(ReflectionType::diffuse, Color3F(0.6f, 0.85f, 1.0f), Color3F(), 1.0f, 0.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position *= 5;
			objects[nr]->triangles[i]->Vertices[j].Position.X -= 4;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Back
	objects[nr]->material = Material(ReflectionType::diffuse, Color3F(0.8f, 1.0f, 0.6f), Color3F(), 1.0f, 0.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position *= 5;
			objects[nr]->triangles[i]->Vertices[j].Position.Z -= 4;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Top
	objects[nr]->material = Material(ReflectionType::diffuse, Color3F(1.0f, 1.0f, 1.0f), Color3F(), 1.0f, 0.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position *= 5;
			objects[nr]->triangles[i]->Vertices[j].Position.Y += 4;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Bottom
	objects[nr]->material = Material(ReflectionType::diffuse, Color3F(1.0f, 1.0f, 1.0f), Color3F(), 1.0f, 0.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position *= 5;
			objects[nr]->triangles[i]->Vertices[j].Position.Y -= 4;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Light
	objects[nr]->material = Material(ReflectionType::diffuse, Color3F(1.0f, 1.0f, 1.0f), Color3F(3.0f, 3.0f, 3.0f), 1.0f, 0.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position /= 1;
			objects[nr]->triangles[i]->Vertices[j].Position.Y += 1.5f;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	lights.push_back(objects[nr]);
	nLight += nTriangles;
	++nr;

	camera = Camera(Vector3F(0.75f, 0, 4.5f), Vector3F::Normalise(Vector3F(-1.0f, 0, -4.5f)), Vector3F(0, 1, 0));

	// hack for initial camera transform bug
	camera.RotateX(0.1f);
	camera.RotateX(-0.1f);
}

void Scene::LoadDefaultScene2()
{
	ObjReader reader;
	Clear();
	unsigned int nTriangles=0, nr=0;

	// Parse central light
	objects = reader.parseFile("../models/cube.obj");

	// Parse left sphere
	objects = reader.parseFile("../models/sphere.obj");

	// Parse right sphere
	objects = reader.parseFile("../models/sphere.obj");

	// Parse front sphere
	objects = reader.parseFile("../models/sphere.obj");

	// Parse back sphere
	objects = reader.parseFile("../models/sphere.obj");

	// Central light
	objects[nr]->material = Material(ReflectionType::refractive, Color3F(1.0f, 1.0f, 1.0f), Color3F(2.5f, 2.5f, 2.5f), 1.5f, 1000000.0f, 1.0f);
	nTriangles = objects[nr]->triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
		objects[nr]->triangles[i]->PreCalc();
	lights.push_back(objects[nr]);
	nLight += nTriangles;
	++nr;

	// Left sphere
	objects[nr]->material = Material(ReflectionType::diffuse, Color3F(1.0f, 0.0f, 0.0f), Color3F(), 1.0f, 1.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position /= 3;
			objects[nr]->triangles[i]->Vertices[j].Position.X -= 1.5f;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Right sphere
	objects[nr]->material = Material(ReflectionType::glossy, Color3F(0.5f, 0.5f, 0.5f), Color3F(), 1.0f, 500.0f, 1.0f);
	nTriangles = objects[nr]->triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position /= 3;
			objects[nr]->triangles[i]->Vertices[j].Position.X += 1.5f;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Front sphere
	objects[nr]->material = Material(ReflectionType::specular, Color3F(1.0f, 1.0f, 1.0f), Color3F(), 1.0f, 1.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position /= 3;
			objects[nr]->triangles[i]->Vertices[j].Position.Z += 1.5f;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Back sphere
	objects[nr]->material = Material(ReflectionType::diffuse, Color3F(0.0f, 0.0f, 1.0f), Color3F(), 1.0f, 1.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position /= 3;
			objects[nr]->triangles[i]->Vertices[j].Position.Z -= 1.5f;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	camera = Camera(Vector3F(-2.5f, 0.75f, 2.5f), Vector3F::Normalise(Vector3F(1, -0.3f, -1)), Vector3F(0, 1, 0));

	// hack for initial camera transform bug
	camera.RotateX(0.1f);
	camera.RotateX(-0.1f);
}

void Scene::LoadDefaultScene3()
{
	ObjReader reader;
	Clear();
	unsigned int nTriangles, nr = 0;

	// Parse left box
	reader.parseFile("../models/cube.obj");

	// Parse right box
	reader.parseFile("../models/cube.obj");

	// Parse right
	reader.parseFile("../models/cube.obj");

	// Parse left
	reader.parseFile("../models/cube.obj");

	// Parse back
	reader.parseFile("../models/cube.obj");

	// Parse top
	reader.parseFile("../models/cube.obj");

	// Parse bottom
	objects = reader.parseFile("../models/cube.obj");

	// Left box
	objects[nr]->material = Material(ReflectionType::diffuse, Color3F(1.0f, 0.0f, 0.0f), Color3F(), 1.0f, 100.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();

	Matrix3x3F mat = Matrix3x3F::CreateRotationY(M_PI_2) * Matrix3x3F::CreateRotationX(M_PI_4);
	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position.X /= 1.25f;
			objects[nr]->triangles[i]->Vertices[j].Position.Z /= 2;

			objects[nr]->triangles[i]->Vertices[j].Position *= mat;
			objects[nr]->triangles[i]->Vertices[j].Normal *= mat;

			objects[nr]->triangles[i]->Vertices[j].Position.X -= 0.75f;
			objects[nr]->triangles[i]->Vertices[j].Position.Y -= 0.97f;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Right box
	objects[nr]->material = Material(ReflectionType::diffuse, Color3F(0.0f, 0.0f, 1.0f), Color3F(), 1.0f, 0.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();

	mat = Matrix3x3F::CreateRotationY(-0.3f);
	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position.X /= 1.25f;
			objects[nr]->triangles[i]->Vertices[j].Position.Y /= 1.5f;
			objects[nr]->triangles[i]->Vertices[j].Position.Z /= 2;

			objects[nr]->triangles[i]->Vertices[j].Position *= mat;
			objects[nr]->triangles[i]->Vertices[j].Normal *= mat;

			objects[nr]->triangles[i]->Vertices[j].Position.X += 0.2f;
			objects[nr]->triangles[i]->Vertices[j].Position.Y -= 7.0f / 6.0f;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Right
	objects[nr]->material = Material(ReflectionType::diffuse, Color3F(1.0f, 1.0f, 1.0f), Color3F(), 1.0f, 0.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position *= 4;
			objects[nr]->triangles[i]->Vertices[j].Position.X += 3.5f;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Left
	objects[nr]->material = Material(ReflectionType::diffuse, Color3F(1.0f, 1.0f, 1.0f), Color3F(), 1.0f, 0.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position *= 4;
			objects[nr]->triangles[i]->Vertices[j].Position.X -= 3.5f;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Back
	objects[nr]->material = Material(ReflectionType::diffuse, Color3F(0.5f, 0.5f, 0.5f), Color3F(), 1.0f, 0.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position *= 4;
			objects[nr]->triangles[i]->Vertices[j].Position.Z -= 3.5f;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Top
	objects[nr]->material = Material(ReflectionType::diffuse, Color3F(1.0f, 1.0f, 1.0f), Color3F(2.0f, 2.0f, 2.0f), 1.0f, 0.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position *= 4;
			objects[nr]->triangles[i]->Vertices[j].Position.Y += 3.5f;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	lights.push_back(objects[nr]);
	nLight += nTriangles;
	++nr;

	// Bottom
	objects[nr]->material = Material(ReflectionType::diffuse, Color3F(0.75f, 0.75f, 0.75f), Color3F(), 1.0f, 0.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position *= 4;
			objects[nr]->triangles[i]->Vertices[j].Position.Y -= 3.5f;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	camera = Camera(Vector3F(0, 0.5f, 3.0f), Vector3F::Normalise(Vector3F(0, -0.5f, -1)), Vector3F(0, 1, 0));

	// hack for initial camera transform bug
	camera.RotateX(0.1f);
	camera.RotateX(-0.1f);
}

void Scene::LoadDefaultScene4()
{
	ObjReader reader;
	Clear();
	unsigned int nTriangles, nr = 0;

	// Parse floor
	reader.parseFile("../models/cube.obj");

	// Parse teapots
	reader.parseFile("../models/teapot.obj");
	reader.parseFile("../models/teapot.obj");
	reader.parseFile("../models/teapot.obj");
	reader.parseFile("../models/teapot.obj");
	reader.parseFile("../models/teapot.obj");

	// Parse light
	objects = reader.parseFile("../models/cube.obj");

	Matrix3x3F mat = Matrix3x3F::CreateRotationY(M_PI + M_PI_4);

	// Floor
	objects[nr]->material = Material(ReflectionType::diffuse, Color3F(0.5f, 0.5f, 0.5f), Color3F(), 1.0f, 0.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position *= 10;
			objects[nr]->triangles[i]->Vertices[j].Position.Y -= 5;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Teapot
	objects[nr]->material = Material(ReflectionType::glossy, Color3F(0.25f, 0.25f, 0.25f), Color3F(), 1.0f, 1.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position *= mat;
			objects[nr]->triangles[i]->Vertices[j].Normal *= mat;

			objects[nr]->triangles[i]->Vertices[j].Position /= 1;
			objects[nr]->triangles[i]->Vertices[j].Position.X -= 2.5f;
			objects[nr]->triangles[i]->Vertices[j].Position.Z -= 1;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Teapot
	objects[nr]->material = Material(ReflectionType::glossy, Color3F(0.25f, 0.25f, 0.25f), Color3F(), 1.0f, 2.5f, 0.0f);
	nTriangles = objects[nr]->triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position *= mat;
			objects[nr]->triangles[i]->Vertices[j].Normal *= mat;

			objects[nr]->triangles[i]->Vertices[j].Position /= 1;
			objects[nr]->triangles[i]->Vertices[j].Position.X -= 1.25f;
			objects[nr]->triangles[i]->Vertices[j].Position.Z -= 0.5f;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Teapot
	objects[nr]->material = Material(ReflectionType::glossy, Color3F(0.25f, 0.25f, 0.25f), Color3F(), 1.0f, 5.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position *= mat;
			objects[nr]->triangles[i]->Vertices[j].Normal *= mat;

			objects[nr]->triangles[i]->Vertices[j].Position /= 1;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Teapot
	objects[nr]->material = Material(ReflectionType::glossy, Color3F(0.25f, 0.25f, 0.25f), Color3F(), 1.0f, 10.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position *= mat;
			objects[nr]->triangles[i]->Vertices[j].Normal *= mat;

			objects[nr]->triangles[i]->Vertices[j].Position /= 1;
			objects[nr]->triangles[i]->Vertices[j].Position.X += 1.25f;
			objects[nr]->triangles[i]->Vertices[j].Position.Z += 0.5f;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Teapot
	objects[nr]->material = Material(ReflectionType::glossy, Color3F(0.25f, 0.25f, 0.25f), Color3F(), 1.0f, 100.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position *= mat;
			objects[nr]->triangles[i]->Vertices[j].Normal *= mat;

			objects[nr]->triangles[i]->Vertices[j].Position /= 1;
			objects[nr]->triangles[i]->Vertices[j].Position.X += 2.5f;
			objects[nr]->triangles[i]->Vertices[j].Position.Z += 1;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Light
	objects[nr]->material = Material(ReflectionType::diffuse, Color3F(1.0f, 1.0f, 1.0f), Color3F(6.0f, 6.0f, 6.0f), 1.0f, 0.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
	{
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position.X *= 10;
			objects[nr]->triangles[i]->Vertices[j].Position.Y += 4;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
		objects[nr]->triangles[i]->PreCalc();
	}
	lights.push_back(objects[nr]);
	nLight += nTriangles;
	++nr;

	camera = Camera(Vector3F(1.34f, 2.405f, 5.25f), Vector3F::Normalise(Vector3F(-0.8f, -2.1f, -5.0f)), Vector3F(0, 1, 0));

	// hack for initial camera transform bug
	camera.RotateX(0.1f);
	camera.RotateX(-0.1f);
}

void Scene::LoadDefaultScene5()
{
	ObjReader reader;
	Clear();
	unsigned int nTriangles, nr = 0;

	// Parse floor
	reader.parseFile("../models/cube.obj");

	// Parse diamonds
	reader.parseFile("../models/diamond.obj");
	reader.parseFile("../models/diamond.obj");
	reader.parseFile("../models/diamond.obj");
	reader.parseFile("../models/diamond.obj");
	reader.parseFile("../models/diamond.obj");

	// Parse light
	objects = reader.parseFile("../models/cube.obj");

	// Floor
	objects[nr]->material = Material(ReflectionType::diffuse, Color3F(0.5f, 0.5f, 0.5f), Color3F(), 1.0f, 0.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position *= 10;
			objects[nr]->triangles[i]->Vertices[j].Position.Y -= 5;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Diamond
	objects[nr]->material = Material(ReflectionType::refractive, Color3F(1.0f, 1.0f, 1.0f), Color3F(), 2.42f, 0.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();
	Matrix3x3F mat = Matrix3x3F::CreateRotationY(M_PI + M_PI_4);
	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position *= mat;
			objects[nr]->triangles[i]->Vertices[j].Normal *= mat;

			objects[nr]->triangles[i]->Vertices[j].Position /= 10;
			objects[nr]->triangles[i]->Vertices[j].Position.X -= 0.25f;
			objects[nr]->triangles[i]->Vertices[j].Position.Z -= 0.25f;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Diamond
	objects[nr]->material = Material(ReflectionType::refractive, Color3F(1.0f, 1.0f, 1.0f), Color3F(), 2.42f, 0.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();
	mat = Matrix3x3F::CreateRotationY(M_PI);
	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position *= mat;
			objects[nr]->triangles[i]->Vertices[j].Normal *= mat;

			objects[nr]->triangles[i]->Vertices[j].Position /= 15;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Diamond
	objects[nr]->material = Material(ReflectionType::refractive, Color3F(1.0f, 1.0f, 1.0f), Color3F(), 2.42f, 0.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();
	mat = Matrix3x3F::CreateRotationY(M_PI_4);
	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position *= mat;
			objects[nr]->triangles[i]->Vertices[j].Normal *= mat;

			objects[nr]->triangles[i]->Vertices[j].Position /= 20;
			objects[nr]->triangles[i]->Vertices[j].Position.X += 0.2f;
			objects[nr]->triangles[i]->Vertices[j].Position.Z -= 0.3f;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Diamond
	objects[nr]->material = Material(ReflectionType::refractive, Color3F(1.0f, 1.0f, 1.0f), Color3F(), 2.42f, 0.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();
	mat = Matrix3x3F::CreateRotationY(0.5f);
	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position *= mat;
			objects[nr]->triangles[i]->Vertices[j].Normal *= mat;

			objects[nr]->triangles[i]->Vertices[j].Position /= 12;
			objects[nr]->triangles[i]->Vertices[j].Position.X += 0.2f;
			objects[nr]->triangles[i]->Vertices[j].Position.Z += 0.3f;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Diamond
	objects[nr]->material = Material(ReflectionType::refractive, Color3F(1.0f, 1.0f, 1.0f), Color3F(), 2.42f, 0.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();
	mat = Matrix3x3F::CreateRotationY(2.0f);
	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position *= mat;
			objects[nr]->triangles[i]->Vertices[j].Normal *= mat;

			objects[nr]->triangles[i]->Vertices[j].Position /= 17;
			objects[nr]->triangles[i]->Vertices[j].Position.X -= 0.15f;
			objects[nr]->triangles[i]->Vertices[j].Position.Z += 0.15f;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Light
	objects[nr]->material = Material(ReflectionType::diffuse, Color3F(1.0f, 1.0f, 1.0f), Color3F(10.0f, 10.0f, 10.0f), 1.0f, 0.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position.X *= 10;
			objects[nr]->triangles[i]->Vertices[j].Position.Y += 4;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	lights.push_back(objects[nr]);
	nLight += nTriangles;
	++nr;

	camera = Camera(Vector3F(1.22f, 0.73f, -0.02f), Vector3F::Normalise(Vector3F(-1.22f, -0.63f, 0.02f)), Vector3F(0, 1, 0));

	// hack for initial camera transform bug
	camera.RotateX(0.1f);
	camera.RotateX(-0.1f);
}

void Scene::LoadDefaultScene6()
{
	ObjReader reader;
	Clear();
	unsigned int nTriangles, nr = 0;

	// Parse floor
	reader.parseFile("../models/cube.obj");

	// Parse diamonds
	reader.parseFile("../models/dragon.obj");
	reader.parseFile("../models/cube.obj");
	reader.parseFile("../models/sphere.obj");
	reader.parseFile("../models/sphere.obj");

	// Parse light
	objects = reader.parseFile("../models/cube.obj");

	// Floor
	objects[nr]->material = Material(ReflectionType::diffuse, Color3F(0.7f, 0.7f, 0.7f), Color3F(), 1.0f, 0.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();
	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position *= 10 * 5;
			objects[nr]->triangles[i]->Vertices[j].Position.Y -= 5 * 5;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Dragon
	objects[nr]->material = Material(ReflectionType::glossy, Color3F((unsigned char)112, 142, 108), Color3F(), 1.0f, 2.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();
	Matrix3x3F mat = Matrix3x3F::CreateRotationY(-M_PI * (5.0f / 16.0f));
	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position *= mat;
			objects[nr]->triangles[i]->Vertices[j].Normal *= mat;

			objects[nr]->triangles[i]->Vertices[j].Position *= 5;
			objects[nr]->triangles[i]->Vertices[j].Position.Y += 0.27f * 5;
			objects[nr]->triangles[i]->Vertices[j].Position.Z -= 0.1f * 5;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Mirror
	objects[nr]->material = Material(ReflectionType::specular, Color3F(0.75f, 0.75f, 0.75f), Color3F(), 1.0f, 0.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();
	mat = Matrix3x3F::CreateRotationY(M_PI_4 * (3.0f / 4.0f));
	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position.X *= 10 * 5;
			objects[nr]->triangles[i]->Vertices[j].Position.Y *= 10 * 5;
			objects[nr]->triangles[i]->Vertices[j].Position.Z *= 0.625f;

			objects[nr]->triangles[i]->Vertices[j].Position *= mat;
			objects[nr]->triangles[i]->Vertices[j].Normal *= mat;

			objects[nr]->triangles[i]->Vertices[j].Position.Y += 5.0f * 5;
			objects[nr]->triangles[i]->Vertices[j].Position.Z += 1.0f * 5;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Glass ball
	objects[nr]->material = Material(ReflectionType::refractive, Color3F(0.7f, 0.7f, 0.7f), Color3F(), 2.5f, 0.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();
	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position /= 0.8f;
			objects[nr]->triangles[i]->Vertices[j].Position.X += 0.3f * 5;
			objects[nr]->triangles[i]->Vertices[j].Position.Y += 0.25f * 5;
			objects[nr]->triangles[i]->Vertices[j].Position.Z += 0.2f * 5;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Shiny ball
	objects[nr]->material = Material(ReflectionType::glossy, Color3F(0.75f, 0.75f, 0.75f), Color3F(), 1.0f, 100.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();
	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position /= 1.0f;
			objects[nr]->triangles[i]->Vertices[j].Position.X += 0.6f * 5;
			objects[nr]->triangles[i]->Vertices[j].Position.Y += 0.2f * 5;
			objects[nr]->triangles[i]->Vertices[j].Position.Z -= 0.1f * 5;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	++nr;

	// Light
	objects[nr]->material = Material(ReflectionType::diffuse, Color3F(1.0f, 1.0f, 1.0f), Color3F(7.0f, 7.0f, 7.0f), 1.0f, 0.0f, 0.0f);
	nTriangles = objects[nr]->triangles.size();
	mat = Matrix3x3F::CreateRotationY(M_PI_4 * (3.0f / 4.0f));
	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr]->triangles[i]->Vertices[j].Position.X *= 5 * 5;
			objects[nr]->triangles[i]->Vertices[j].Position.Y /= 0.4f;
			objects[nr]->triangles[i]->Vertices[j].Position.Z *= 2 * 5;

			objects[nr]->triangles[i]->Vertices[j].Position *= mat;
			objects[nr]->triangles[i]->Vertices[j].Normal *= mat;

			objects[nr]->triangles[i]->Vertices[j].Position.X -= 0.5f * 5;
			objects[nr]->triangles[i]->Vertices[j].Position.Y += 4.0f * 5;
			objects[nr]->triangles[i]->Vertices[j].Position.Z -= 0.5f * 5;
		}
		objects[nr]->triangles[i]->PreCalc();
	}
	lights.push_back(objects[nr]);
	nLight += nTriangles;
	++nr;

	camera = Camera(Vector3F(0.3f * 5, 1.25f * 5, -1.5f * 5), Vector3F::Normalise(Vector3F(0, -1.2f, 1.75f)), Vector3F(0, 1, 0));

	// hack for initial camera transform bug
	camera.RotateX(0.1f);
	camera.RotateX(-0.1f);
}

void Scene::Clear()
{
	for (unsigned int i = 0; i < objects.size(); i++)
		delete objects[i];
	objects.clear();
	lights.clear();
	nLight = 0;
	nLight8 = 0;
}