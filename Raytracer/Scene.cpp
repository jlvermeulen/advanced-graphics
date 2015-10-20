#define _USE_MATH_DEFINES
#define RUSSIAN_ROULETTE_PROBABILITY 0.90
#define MIN_PATH_LENGTH 2
#define TREE 1 // 1 = bhv 0 = octtree

#include "Scene.h"

#include "Intersections.h"
#include "ObjReader.h"
#include <algorithm>
#include <math.h>
#include <random>
#include "IntersectionInfo.h"

Scene::Scene() :
  gen(std::mt19937(std::random_device()())),
  dist(std::uniform_real_distribution<double>(0, 1))
{
	// checkerboard
	//Vertex v1(Vector3D(-100, -0.5, -100), Vector3D(0, 1, 0), ColorD(1.0, 1.0, 1.0), Vector3D());
	//Vertex v2(Vector3D(100, -0.5, -100), Vector3D(0, 1, 0), ColorD(1.0, 1.0, 1.0), Vector3D());
	//Vertex v3(Vector3D(-100, -0.5, 100), Vector3D(0, 1, 0), ColorD(1.0, 1.0, 1.0), Vector3D());
	//Vertex v4(Vector3D(100, -0.5, 100), Vector3D(0, 1, 0), ColorD(1.0, 1.0, 1.0), Vector3D());

	//std::deque<Triangle> tris;
	//tris.push_back(Triangle(v1, v3, v2));
	//tris.push_back(Triangle(v4, v2, v3));

	//checkerboard = Object(tris, Material(ReflectionType::diffuse, ColorD(1.0, 1.0, 1.0), ColorD(), 1, 0));

	// Add lights
	//lights.push_back(Light(Vector3D(-3.0, -5.0, -4.0), ColorD(25.0, 25.0, 25.0)));
	//lights.push_back(Light(Vector3D(3.0, 5.0, 4.0), ColorD(25.0, 25.0, 25.0)));

	//camera = Camera(Vector3D(0, 15.0, -10.0), Vector3D::Normalise(Vector3D(0, -0.5, 1)), Vector3D(0, 1, 0));
}

Scene::~Scene()
{

}

bool Scene::Render(uchar* imageData, int minTriangles, int maxDepth, int samplesPerPixel, double sigma, bool useDoF)
{
	// Instantiate octrees
	for (Object& obj : objects)
	{
#if TREE == 0
		obj.ConstructOctree(minTriangles, maxDepth);
#else
	obj.ConstructBVH(minTriangles);
#endif
}
	std::pair<ColorD, double>* samples = new std::pair<ColorD, double>[camera.Width * camera.Height];
	TracePixels(samples, samplesPerPixel, sigma, useDoF);

	#pragma omp parallel for
	for (int x = 0; x < camera.Width; ++x)
	{
		for (int y = 0; y < camera.Height; ++y)
		{
			std::pair<ColorD, double> sample = samples[y * camera.Width + x];
			ColorD color = sample.first / sample.second;
			color.Clip();

			int offset = (y * camera.Width + x) * 4;

			// For each color channel in reversed order (i.e. blue-green-red-alpha)
			imageData[offset] = (uchar) (color.B * 255.0);
			imageData[offset + 1] = (uchar) (color.G * 255.0);
			imageData[offset + 2] = (uchar) (color.R * 255.0);
			imageData[offset + 3] = 255;
		}
	}

	delete [] samples;

	return true;
}

//--------------------------------------------------------------------------------
void Scene::TracePixels(std::pair<ColorD, double>* pixelData, int samplesPerPixel, double sigma, bool useDoF)
{
	double tanHalfFovY = tan(camera.FovY() / 360 * M_PI);
	double tanHalfFovX = tanHalfFovY * camera.Width / camera.Height;

	double left = -tanHalfFovX;
	double right = tanHalfFovX;
	double top = tanHalfFovY;
	double bottom = -tanHalfFovY;

	// Calculate pixel rays
	for (int o = 0; o < 3; ++o) // prevent concurrent access of same array indices
		#pragma omp parallel for schedule(dynamic)
		for (int x = o; x < camera.Width; x += 3)
			for (int y = 0; y < camera.Height; ++y)
				for (int r = 0; r < samplesPerPixel; ++r)
				{
					double rX = dist(gen);
					double rY = dist(gen);

					double a = left + (right - left) * (x + rX) / camera.Width;
					double b = top + (bottom - top) * (y + rY) / camera.Height;

					Vector3D origin = camera.Eye();
					Vector3D direction = Vector3D::Normalise(camera.Focus() + a * camera.Right() + b * camera.Up());

					if (useDoF)
					{
						Vector3D focalPoint = origin + direction * camera.FocalDepth;

						// Get uniformly distributed square [-1,1] x [-1,1]
						float angle = dist(gen) * 2 * M_PI;
						float radius = dist(gen);

						double aX = cos(angle) * radius * camera.Aperture;
						double aY = sin(angle) * radius * camera.Aperture;

						Vector3D apertureOffset(aX, aY, 0.0);

						origin += apertureOffset;

						direction *= camera.FocalDepth;
						direction -= apertureOffset;
						direction.Normalise();

						//double uX = 2.0 * dist(gen) - 1.0;
						//double uY = 2.0 * dist(gen) - 1.0;

						// Get uniformly distributed circle with lens radius
						//double lX = camera.Aperture * uX * sqrt(1 - uY * uY / 2);
						//double lY = camera.Aperture * uY * sqrt(1 - uX * uX / 2);

						//origin += lX * camera.Right() + lY * camera.Up();

						//direction = Vector3D::Normalise(focalPoint - origin);
					}

					Ray cameraRay(origin, direction);

					//for (unsigned int c = 0; c < 3; ++c)
					{
						ColorD value = TraceRay(cameraRay);

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

								double weight = GaussianWeight(i - rX + 0.5, j - rY + 0.5, sigma);

								int index = yy * camera.Width + xx;
								std::pair<ColorD, double>& data = pixelData[index];

								data.first += value * weight;
								data.second += weight;
							}
						}
					}
				}
}

//--------------------------------------------------------------------------------
double Scene::GaussianWeight(double dx, double dy, double sigma) const
{
	return exp(-(dx * dx + dy * dy) / (2 * sigma * sigma));
}

//--------------------------------------------------------------------------------
ColorD Scene::TraceRay(const Ray& ray)
{
	Material hitMaterial;
	Triangle hitTriangle;
	double hitTime;

	if (!FirstHitInfo(ray, hitTime, hitTriangle, hitMaterial))
		return ColorD();

	return ComputeRadiance(ray.Origin + hitTime * ray.Direction, ray.Direction, hitTriangle, hitMaterial, 0);
}

//--------------------------------------------------------------------------------
ColorD Scene::ComputeRadiance(const Vector3D& point, const Vector3D& in, const Triangle& triangle, const Material& material, unsigned int depth)
{
	const Vector3D& normal = triangle.surfaceNormal(point);
	ColorD c = material.emission * pow(abs(Vector3D::Dot(in, normal)), 1.0 / 10.0);
	c += DirectIllumination(point, normal, material);
	c += IndirectIllumination(point, in, normal, material, depth);
	return c;
}

ColorD Scene::DirectIllumination(const Vector3D& point, const Vector3D& normal, const Material& material)
{
	std::pair<Ray, Triangle> sample = SampleLight(point);

	Material hitMaterial;
	Triangle hitTriangle;
	double hitTime;

	Ray ray(sample.first.Origin/* + 0.005 * normal*/, sample.first.Direction);
	if (!FirstHitInfo(ray, hitTime, hitTriangle, hitMaterial) || (hitTriangle != sample.second && hitMaterial.reflType != ReflectionType::refractive))
		return ColorD();

	if (hitMaterial.reflType == ReflectionType::refractive && !hitMaterial.emission.IsSignificant())
	{
		Vector3D hitPoint =  ray.Origin + ray.Direction * hitTime;
		return DirectIllumination(hitPoint, hitTriangle.surfaceNormal(hitPoint), material);
	}

	return hitMaterial.emission * material.color;
}

ColorD Scene::IndirectIllumination(Vector3D point, const Vector3D& in, const Vector3D& normal, const Material& material, unsigned int depth)
{
	if (depth > MIN_PATH_LENGTH && dist(gen) > RUSSIAN_ROULETTE_PROBABILITY)
		return ColorD();

	Ray ray(point, in);

	double c = 1;// / (2 * M_PI);
	ColorD value(c, c, c);

	// Sample unit sphere
	double phi = 2.0 * M_PI * dist(gen);

	double cosTheta = pow(1.0 - dist(gen), 1.0 / (1.0 + material.specularExponent));
	double sinTheta = sqrt(1.0 - cosTheta * cosTheta);

	double x = sinTheta * cos(phi);
	double y = sinTheta * sin(phi);
	double z = cosTheta;

	Vector3D hemi(x, y, z);

	//double z = dist(gen) * 2 - 1, theta = dist(gen) * M_PI * 2, r = sqrt(1 - z * z);

	//Vector3D hemi(r * cos(theta), r * sin(theta), z);
	//if (hemi.Dot(normal) < 0) // flip if "behind" normal
	//  hemi *= -1;

	if (material.reflType == ReflectionType::specular)
	{
		if (Vector3D::Dot(ray.Direction, normal) < 0)
			ray.Reflect(point, normal);
		//value = ColorD(1.0, 1.0, 1.0);
	}
	else if (material.reflType == ReflectionType::glossy)
	{
		ray.Reflect(point, normal);

		Vector3D w = ray.Direction;

		Vector3D u = Vector3D::Cross(Vector3D(0.00424, 1, 0.00764), w);
		u.Normalise();

		Vector3D v = Vector3D::Cross(u, w);

		Vector3D reflDir = hemi.X * u + hemi.Y * v + hemi.Z * w;

		if (Vector3D::Dot(reflDir, normal) < 0)
			reflDir *= -1;

		ray = Ray(point, reflDir);
		value *= material.color * Vector3D::Dot(normal, w);
	}
	else if (material.reflType == ReflectionType::diffuse)
	{
		double d = Vector3D::Dot(hemi, normal);
		if (d < 0) // flip if "behind" normal
		{
			hemi *= -1;
			d *= -1;
		}

		ray = Ray(point, hemi);
		point += 0.005 * normal;
		value *= material.color * d;
	}
	else if (material.reflType == ReflectionType::refractive)
		ray.Refract(point, normal, 1.0, material.refrIndex);

	Material hitMaterial;
	Triangle hitTriangle;
	double hitTime;

	if (!FirstHitInfo(ray, hitTime, hitTriangle, hitMaterial))
		return ColorD();

	Vector3D hitPoint = ray.Origin + hitTime * ray.Direction;
	return value * ComputeRadiance(hitPoint, ray.Direction, hitTriangle, hitMaterial, depth + 1) / (depth > MIN_PATH_LENGTH ? RUSSIAN_ROULETTE_PROBABILITY : 1.0);
}

//--------------------------------------------------------------------------------
bool Scene::FirstHitInfo(const Ray& ray, double& time, Triangle& triangle, Material& mat) const
{
	
	time = std::numeric_limits<double>::max();
	bool hit = false;
#if TREE==0
	for (const Object& obj : objects)
	{
		Triangle tri;
		double t = std::numeric_limits<double>::max();

		if (obj.octree->Query(ray, tri, t) && t < time)
		{
			mat = obj.material;
			triangle = tri;
			time = t;
			hit = true;
		}
	}
#else
	for (const Object& obj : objects)
	{
		Triangle tri;
		//double t = std::numeric_limits<double>::max();
		IntersectionInfo io;
		if (obj.bvh->getIntersection(ray,&io,false) && io.t < time)
		{
			mat = obj.material;
			triangle = *io.triangle;
			time = io.t;
			hit = true;
		}		
	}
#endif
	return hit;
}

//--------------------------------------------------------------------------------
std::pair<Ray, Triangle> Scene::SampleLight(const Vector3D& hitPoint)
{
	std::deque<std::pair<double, const Triangle&>> flux;
	double totalflux = 0;
	// loop over all lightsources
	for (const Object& o : lights)
		for (const Triangle& t : o.triangles)
		{
			Vector3D outgoingRay = hitPoint - t.Center;
			double distance = outgoingRay.LengthSquared();
			Vector3D normal = t.surfaceNormal(t.Center);
			outgoingRay.Normalise();

			double f = t.Area * o.material.emission.R * std::max(0.0, Vector3D::Dot(outgoingRay, normal)) / distance;
			flux.push_back(std::pair<double, const Triangle&>(f, t));
			totalflux += f;
		}

	// divide all flux values by the total flux gives the probabilities
	for (unsigned int l = 0; l < flux.size(); l++)
		flux[l].first /= totalflux;

	// get the light using the probabilities and a random float between 0 and 1
	double r = dist(gen);
	std::pair<double, Triangle> pair; int x = -1;
	for (std::pair<double, const Triangle&> p : flux)
	{
		if (p.first >= r)
		{
			pair = p;
			x = 0;
			break;
		}
		else
			r -= p.first;
	}

	// create random point within the lightsource and return ray towards it.
	double u = dist(gen), v = dist(gen);
	if (u + v > 1)
	{
		u = 1 - u;
		v = 1 - v;
	}

	Vector3D v1 = pair.second.Vertices[1].Position - pair.second.Vertices[0].Position;
	Vector3D v2 = pair.second.Vertices[2].Position - pair.second.Vertices[0].Position;

	return std::pair<Ray, Triangle>(Ray(hitPoint, Vector3D::Normalise(v1 * u + v2 * v + pair.second.Vertices[0].Position - hitPoint)), pair.second);
}

//--------------------------------------------------------------------------------
void Scene::LoadDefaultScene()
{
	ObjReader reader;
	objects.clear();
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
	objects[nr].material = Material(ReflectionType::glossy, ColorD(0.75, 0.25, 0.25), ColorD(), 1.0, 100.0, 0.0);
	nTriangles = objects[nr].triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr].triangles[i].Vertices[j].Position /= 2;

			objects[nr].triangles[i].Vertices[j].Position.X += 0.8;
		}
	++nr;

	// Left sphere
	objects[nr].material = Material(ReflectionType::specular, ColorD(0.1, 0.1, 0.1), ColorD(), 1.0, 0.0, 0.0);
	nTriangles = objects[nr].triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr].triangles[i].Vertices[j].Position /= 3.5;

			objects[nr].triangles[i].Vertices[j].Position.X += 0.3;
			objects[nr].triangles[i].Vertices[j].Position.Y += 0.8;
			objects[nr].triangles[i].Vertices[j].Position.Z -= 1.0;
		}
	++nr;

	// Left sphere
	objects[nr].material = Material(ReflectionType::refractive, ColorD(0.0, 0.0, 0.0), ColorD(), 1.5, 0.0, 0.0);
	nTriangles = objects[nr].triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr].triangles[i].Vertices[j].Position /= 1.5;

			objects[nr].triangles[i].Vertices[j].Position.X -= 0.7;
			objects[nr].triangles[i].Vertices[j].Position.Y -= 5.0 / 6.0;
			objects[nr].triangles[i].Vertices[j].Position.Z -= 0.7;
		}
	++nr;

	// Left sphere
	objects[nr].material = Material(ReflectionType::glossy, ColorD(0.25, 0.25, 0.75), ColorD(), 1.0, 10.0, 0.0);
	nTriangles = objects[nr].triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr].triangles[i].Vertices[j].Position /= 3;

			objects[nr].triangles[i].Vertices[j].Position.X -= 0.75;
			objects[nr].triangles[i].Vertices[j].Position.Y += 0.75;
			objects[nr].triangles[i].Vertices[j].Position.Z += 2.0;
		}
	++nr;

	// Left sphere
	objects[nr].material = Material(ReflectionType::glossy, ColorD(0.5, 0.0, 0.0), ColorD(), 1.0, 0.0, 0.0);
	nTriangles = objects[nr].triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr].triangles[i].Vertices[j].Position /= 4;

			objects[nr].triangles[i].Vertices[j].Position.X += 0.5;
			objects[nr].triangles[i].Vertices[j].Position.Y -= 1.25;
			objects[nr].triangles[i].Vertices[j].Position.Z += 0.5;
		}
	++nr;

	// Right
	objects[nr].material = Material(ReflectionType::diffuse, ColorD(1.0, 0.6, 0.6), ColorD(), 1.0, 0.0, 0.0);
	nTriangles = objects[nr].triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr].triangles[i].Vertices[j].Position *= 5;
			objects[nr].triangles[i].Vertices[j].Position.X += 4;
		}
	++nr;

	// Left
	objects[nr].material = Material(ReflectionType::diffuse, ColorD(0.6, 0.85, 1.0), ColorD(), 1.0, 0.0, 0.0);
	nTriangles = objects[nr].triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr].triangles[i].Vertices[j].Position *= 5;
			objects[nr].triangles[i].Vertices[j].Position.X -= 4;
		}
	++nr;

	// Back
	objects[nr].material = Material(ReflectionType::diffuse, ColorD(0.8, 1.0, 0.6), ColorD(), 1.0, 0.0, 0.0);
	nTriangles = objects[nr].triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr].triangles[i].Vertices[j].Position *= 5;
			objects[nr].triangles[i].Vertices[j].Position.Z -= 4;
		}
	++nr;

	// Top
	objects[nr].material = Material(ReflectionType::diffuse, ColorD(1.0, 1.0, 1.0), ColorD(), 1.0, 0.0, 0.0);
	nTriangles = objects[nr].triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr].triangles[i].Vertices[j].Position *= 5;
			objects[nr].triangles[i].Vertices[j].Position.Y += 4;
		}
		objects[nr].triangles[i].CalculateArea();
		objects[nr].triangles[i].CalculateCenter();
	}
	lights.push_back(objects[nr]);
	++nr;

	// Bottom
	objects[nr].material = Material(ReflectionType::diffuse, ColorD(1.0, 1.0, 1.0), ColorD(), 1.0, 0.0, 0.0);
	nTriangles = objects[nr].triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr].triangles[i].Vertices[j].Position *= 5;
			objects[nr].triangles[i].Vertices[j].Position.Y -= 4;
		}
	++nr;

	// Light
	objects[nr].material = Material(ReflectionType::diffuse, ColorD(1.0, 1.0, 1.0), ColorD(10.0, 10.0, 10.0), 1.0, 0.0, 0.0);
	nTriangles = objects[nr].triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr].triangles[i].Vertices[j].Position /= 1;
			objects[nr].triangles[i].Vertices[j].Position.Y += 1.5;
		}
	++nr;

	camera = Camera(Vector3D(0.75, 0, 4.5), Vector3D::Normalise(Vector3D(-1.0, 0, -4.5)), Vector3D(0, 1, 0));

	// hack for initial camera transform bug
	camera.RotateX(0.1f);
	camera.RotateX(-0.1f);
}

void Scene::LoadDefaultScene2()
{
	ObjReader reader;
	objects.clear();
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
	objects[nr].material = Material(ReflectionType::refractive, ColorD(0.0, 0.0, 0.0), ColorD(0.5, 0.5, 0.5), 1.5, 1000000.0, 1.0);
	nTriangles = objects[nr].triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
	{
		objects[nr].triangles[i].CalculateArea();
		objects[nr].triangles[i].CalculateCenter();
	}
	lights.push_back(objects[nr]);
	++nr;

	// Left sphere
	objects[nr].material = Material(ReflectionType::diffuse, ColorD(1.0, 0.0, 0.0), ColorD(), 1.0, 1.0, 0.0);
	nTriangles = objects[nr].triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr].triangles[i].Vertices[j].Position /= 3;
			objects[nr].triangles[i].Vertices[j].Position.X -= 1.5;
		}
	++nr;

	// Right sphere
	objects[nr].material = Material(ReflectionType::glossy, ColorD(0.5, 0.5, 0.5), ColorD(), 1.0, 500.0, 1.0);
	nTriangles = objects[nr].triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr].triangles[i].Vertices[j].Position /= 3;
			objects[nr].triangles[i].Vertices[j].Position.X += 1.5;
		}
	++nr;

	// Front sphere
	objects[nr].material = Material(ReflectionType::specular, ColorD(0.1, 0.1, 0.1), ColorD(), 1.0, 1.0, 0.0);
	nTriangles = objects[nr].triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr].triangles[i].Vertices[j].Position /= 3;
			objects[nr].triangles[i].Vertices[j].Position.Z += 1.5;
		}
	++nr;

	// Back sphere
	objects[nr].material = Material(ReflectionType::diffuse, ColorD(0.0, 0.0, 1.0), ColorD(), 1.0, 1.0, 0.0);
	nTriangles = objects[nr].triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr].triangles[i].Vertices[j].Position /= 3;
			objects[nr].triangles[i].Vertices[j].Position.Z -= 1.5;
		}
	++nr;

	camera = Camera(Vector3D(-2.5, 0.75, 2.5), Vector3D::Normalise(Vector3D(1, -0.3, -1)), Vector3D(0, 1, 0));

	// hack for initial camera transform bug
	camera.RotateX(0.1f);
	camera.RotateX(-0.1f);
}

void Scene::LoadDefaultScene3()
{
	ObjReader reader;
	objects.clear();
	unsigned int nTriangles, nr = 0;

	// Parse left box
	//reader.parseFile("../models/cube.obj");

	// Parse right box
//	reader.parseFile("../models/cube.obj");

	// Parse right
//	reader.parseFile("../models/cube.obj");

	// Parse left
//	reader.parseFile("../models/cube.obj");

	// Parse back
	//reader.parseFile("../models/cube.obj");

	// Parse top
//	reader.parseFile("../models/cube.obj");

	// Parse bottom
	objects = reader.parseFile("../models/cube.obj");

	// Left box
	//objects[nr].material = Material(ReflectionType::glossy, ColorD(1.0, 0.0, 0.0), ColorD(), 1.0, 100.0, 0.0);
	//nTriangles = objects[nr].triangles.size();

	//Matrix3x3D mat = Matrix3x3D::CreateRotationY(M_PI_2) * Matrix3x3D::CreateRotationX(M_PI_4);
	//for (unsigned int i = 0; i < nTriangles; ++i)
	//	for (unsigned int j = 0; j < 3; ++j)
	//	{
	//		objects[nr].triangles[i].Vertices[j].Position.X /= 1.25;
	//		objects[nr].triangles[i].Vertices[j].Position.Z /= 2;

	//		objects[nr].triangles[i].Vertices[j].Position *= mat;
	//		objects[nr].triangles[i].Vertices[j].Normal *= mat;

	//		objects[nr].triangles[i].Vertices[j].Position.X -= 0.75;
	//		objects[nr].triangles[i].Vertices[j].Position.Y -= 0.97;
	//	}
	//++nr;

	//// Right box
	//objects[nr].material = Material(ReflectionType::diffuse, ColorD(0.0, 0.0, 1.0), ColorD(), 1.0, 0.0, 0.0);
	//nTriangles = objects[nr].triangles.size();

	//mat = Matrix3x3D::CreateRotationY(-0.3);
	//for (unsigned int i = 0; i < nTriangles; ++i)
	//	for (unsigned int j = 0; j < 3; ++j)
	//	{
	//		objects[nr].triangles[i].Vertices[j].Position.X /= 1.25;
	//		objects[nr].triangles[i].Vertices[j].Position.Y /= 1.5;
	//		objects[nr].triangles[i].Vertices[j].Position.Z /= 2;

	//		objects[nr].triangles[i].Vertices[j].Position *= mat;
	//		objects[nr].triangles[i].Vertices[j].Normal *= mat;

	//		objects[nr].triangles[i].Vertices[j].Position.X += 0.2;
	//		objects[nr].triangles[i].Vertices[j].Position.Y -= 7.0 / 6.0;
	//	}
	//++nr;

	//// Right
	//objects[nr].material = Material(ReflectionType::diffuse, ColorD(1.0, 1.0, 1.0), ColorD(), 1.0, 0.0, 0.0);
	//nTriangles = objects[nr].triangles.size();

	//for (unsigned int i = 0; i < nTriangles; ++i)
	//	for (unsigned int j = 0; j < 3; ++j)
	//	{
	//		objects[nr].triangles[i].Vertices[j].Position *= 5;
	//		objects[nr].triangles[i].Vertices[j].Position.X += 4;
	//	}
	//++nr;

	//// Left
	//objects[nr].material = Material(ReflectionType::diffuse, ColorD(1.0, 1.0, 1.0), ColorD(), 1.0, 0.0, 0.0);
	//nTriangles = objects[nr].triangles.size();

	//for (unsigned int i = 0; i < nTriangles; ++i)
	//	for (unsigned int j = 0; j < 3; ++j)
	//	{
	//		objects[nr].triangles[i].Vertices[j].Position *= 5;
	//		objects[nr].triangles[i].Vertices[j].Position.X -= 4;
	//	}
	//++nr;

	// Back
	objects[nr].material = Material(ReflectionType::diffuse, ColorD(0.5, 0.5, 0.5), ColorD(), 1.0, 0.0, 0.0);
	nTriangles = objects[nr].triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr].triangles[i].Vertices[j].Position *= 5;
			objects[nr].triangles[i].Vertices[j].Position.Z -= 2;
		}
	++nr;

	//// Top
	//objects[nr].material = Material(ReflectionType::diffuse, ColorD(1.0, 1.0, 1.0), ColorD(1.0, 1.0, 1.0), 1.0, 0.0, 0.0);
	//nTriangles = objects[nr].triangles.size();

	//for (unsigned int i = 0; i < nTriangles; ++i)
	//{
	//	for (unsigned int j = 0; j < 3; ++j)
	//	{
	//		objects[nr].triangles[i].Vertices[j].Position *= 5;
	//		objects[nr].triangles[i].Vertices[j].Position.Y += 4;
	//	}
	//	objects[nr].triangles[i].CalculateArea();
	//	objects[nr].triangles[i].CalculateCenter();
	//}
	//lights.push_back(objects[nr]);
	//++nr;

	//// Bottom
	//objects[nr].material = Material(ReflectionType::diffuse, ColorD(0.75, 0.75, 0.75), ColorD(), 1.0, 0.0, 0.0);
	//nTriangles = objects[nr].triangles.size();

	//for (unsigned int i = 0; i < nTriangles; ++i)
	//	for (unsigned int j = 0; j < 3; ++j)
	//	{
	//		objects[nr].triangles[i].Vertices[j].Position *= 5;
	//		objects[nr].triangles[i].Vertices[j].Position.Y -= 4;
	//	}
	//++nr;

	camera = Camera(Vector3D(0, 0.5, 3.0), Vector3D::Normalise(Vector3D(0, -0.5, -1)), Vector3D(0, 1, 0));

	// hack for initial camera transform bug
	camera.RotateX(0.1f);
	camera.RotateX(-0.1f);
}

void Scene::LoadDefaultScene4()
{
	ObjReader reader;
	objects.clear();
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

	Matrix3x3D mat = Matrix3x3D::CreateRotationY(M_PI + M_PI_4);

	// Left box
	objects[nr].material = Material(ReflectionType::glossy, ColorD(0.5, 0.5, 0.5), ColorD(), 1.0, 0.0, 0.0);
	nTriangles = objects[nr].triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr].triangles[i].Vertices[j].Position *= 10;
			objects[nr].triangles[i].Vertices[j].Position.Y -= 5;
		}
	++nr;

	// Teapot
	objects[nr].material = Material(ReflectionType::glossy, ColorD(0.25, 0.25, 0.25), ColorD(), 1.0, 0.0, 0.0);
	nTriangles = objects[nr].triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr].triangles[i].Vertices[j].Position *= mat;
			objects[nr].triangles[i].Vertices[j].Normal *= mat;

			objects[nr].triangles[i].Vertices[j].Position /= 1;
			objects[nr].triangles[i].Vertices[j].Position.X -= 2.5;
			objects[nr].triangles[i].Vertices[j].Position.Z -= 1;
		}
	++nr;

	// Teapot
	objects[nr].material = Material(ReflectionType::glossy, ColorD(0.25, 0.25, 0.25), ColorD(), 1.0, 2.5, 0.0);
	nTriangles = objects[nr].triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr].triangles[i].Vertices[j].Position *= mat;
			objects[nr].triangles[i].Vertices[j].Normal *= mat;

			objects[nr].triangles[i].Vertices[j].Position /= 1;
			objects[nr].triangles[i].Vertices[j].Position.X -= 1.25;
			objects[nr].triangles[i].Vertices[j].Position.Z -= 0.5;
		}
	++nr;

	// Teapot
	objects[nr].material = Material(ReflectionType::glossy, ColorD(0.25, 0.25, 0.25), ColorD(), 1.0, 5.0, 0.0);
	nTriangles = objects[nr].triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr].triangles[i].Vertices[j].Position *= mat;
			objects[nr].triangles[i].Vertices[j].Normal *= mat;

			objects[nr].triangles[i].Vertices[j].Position /= 1;
		}
	++nr;

	// Teapot
	objects[nr].material = Material(ReflectionType::glossy, ColorD(0.25, 0.25, 0.25), ColorD(), 1.0, 10.0, 0.0);
	nTriangles = objects[nr].triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr].triangles[i].Vertices[j].Position *= mat;
			objects[nr].triangles[i].Vertices[j].Normal *= mat;

			objects[nr].triangles[i].Vertices[j].Position /= 1;
			objects[nr].triangles[i].Vertices[j].Position.X += 1.25;
			objects[nr].triangles[i].Vertices[j].Position.Z += 0.5;
		}
	++nr;

	// Teapot
	objects[nr].material = Material(ReflectionType::glossy, ColorD(0.25, 0.25, 0.25), ColorD(), 1.0, 100.0, 0.0);
	nTriangles = objects[nr].triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr].triangles[i].Vertices[j].Position *= mat;
			objects[nr].triangles[i].Vertices[j].Normal *= mat;

			objects[nr].triangles[i].Vertices[j].Position /= 1;
			objects[nr].triangles[i].Vertices[j].Position.X += 2.5;
			objects[nr].triangles[i].Vertices[j].Position.Z += 1;
		}
	++nr;

	// Light
	objects[nr].material = Material(ReflectionType::diffuse, ColorD(1.0, 1.0, 1.0), ColorD(10.0, 10.0, 10.0), 1.0, 0.0, 0.0);
	nTriangles = objects[nr].triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr].triangles[i].Vertices[j].Position.X *= 10;
			objects[nr].triangles[i].Vertices[j].Position.Y += 4;
		}
	++nr;

	camera = Camera(Vector3D(1.34, 2.405, 5.25), Vector3D::Normalise(Vector3D(-0.8, -2.1, -5.0)), Vector3D(0, 1, 0));

	// hack for initial camera transform bug
	camera.RotateX(0.1f);
	camera.RotateX(-0.1f);
}

void Scene::LoadDefaultScene5()
{
	ObjReader reader;
	objects.clear();
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

	// Left box
	objects[nr].material = Material(ReflectionType::diffuse, ColorD(0.5, 0.5, 0.5), ColorD(), 1.0, 0.0, 0.0);
	nTriangles = objects[nr].triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr].triangles[i].Vertices[j].Position *= 10;
			objects[nr].triangles[i].Vertices[j].Position.Y -= 5;
		}
	++nr;

	// Diamond
	objects[nr].material = Material(ReflectionType::refractive, ColorD(), ColorD(), 2.42, 0.0, 0.0);
	nTriangles = objects[nr].triangles.size();
	Matrix3x3D mat = Matrix3x3D::CreateRotationY(M_PI + M_PI_4);
	for (unsigned int i = 0; i < nTriangles; ++i)
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr].triangles[i].Vertices[j].Position *= mat;
			objects[nr].triangles[i].Vertices[j].Normal *= mat;

			objects[nr].triangles[i].Vertices[j].Position /= 10;
			objects[nr].triangles[i].Vertices[j].Position.X -= 0.25;
			objects[nr].triangles[i].Vertices[j].Position.Z -= 0.25;
		}
	++nr;

	// Diamond
	objects[nr].material = Material(ReflectionType::refractive, ColorD(), ColorD(), 2.42, 0.0, 0.0);
	nTriangles = objects[nr].triangles.size();
	mat = Matrix3x3D::CreateRotationY(M_PI);
	for (unsigned int i = 0; i < nTriangles; ++i)
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr].triangles[i].Vertices[j].Position *= mat;
			objects[nr].triangles[i].Vertices[j].Normal *= mat;

			objects[nr].triangles[i].Vertices[j].Position /= 15;
		}
	++nr;

	// Diamond
	objects[nr].material = Material(ReflectionType::refractive, ColorD(), ColorD(), 2.42, 0.0, 0.0);
	nTriangles = objects[nr].triangles.size();
	mat = Matrix3x3D::CreateRotationY(M_PI_4);
	for (unsigned int i = 0; i < nTriangles; ++i)
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr].triangles[i].Vertices[j].Position *= mat;
			objects[nr].triangles[i].Vertices[j].Normal *= mat;

			objects[nr].triangles[i].Vertices[j].Position /= 20;
			objects[nr].triangles[i].Vertices[j].Position.X += 0.2;
			objects[nr].triangles[i].Vertices[j].Position.Z -= 0.3;
		}
	++nr;

	// Diamond
	objects[nr].material = Material(ReflectionType::refractive, ColorD(), ColorD(), 2.42, 0.0, 0.0);
	nTriangles = objects[nr].triangles.size();
	mat = Matrix3x3D::CreateRotationY(0.5);
	for (unsigned int i = 0; i < nTriangles; ++i)
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr].triangles[i].Vertices[j].Position *= mat;
			objects[nr].triangles[i].Vertices[j].Normal *= mat;

			objects[nr].triangles[i].Vertices[j].Position /= 12;
			objects[nr].triangles[i].Vertices[j].Position.X += 0.2;
			objects[nr].triangles[i].Vertices[j].Position.Z += 0.3;
		}
	++nr;

	// Diamond
	objects[nr].material = Material(ReflectionType::refractive, ColorD(), ColorD(), 2.42, 0.0, 0.0);
	nTriangles = objects[nr].triangles.size();
	mat = Matrix3x3D::CreateRotationY(2.0);
	for (unsigned int i = 0; i < nTriangles; ++i)
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr].triangles[i].Vertices[j].Position *= mat;
			objects[nr].triangles[i].Vertices[j].Normal *= mat;

			objects[nr].triangles[i].Vertices[j].Position /= 17;
			objects[nr].triangles[i].Vertices[j].Position.X -= 0.15;
			objects[nr].triangles[i].Vertices[j].Position.Z += 0.15;
		}
	++nr;

	// Light
	objects[nr].material = Material(ReflectionType::diffuse, ColorD(1.0, 1.0, 1.0), ColorD(10.0, 10.0, 10.0), 1.0, 0.0, 0.0);
	nTriangles = objects[nr].triangles.size();

	for (unsigned int i = 0; i < nTriangles; ++i)
		for (unsigned int j = 0; j < 3; ++j)
		{
			objects[nr].triangles[i].Vertices[j].Position.X *= 10;
			objects[nr].triangles[i].Vertices[j].Position.Y += 4;
		}
	++nr;

	camera = Camera(Vector3D(1.22, 0.73, -0.02), Vector3D::Normalise(Vector3D(-1.22, -0.63, 0.02)), Vector3D(0, 1, 0));

	// hack for initial camera transform bug
	camera.RotateX(0.1f);
	camera.RotateX(-0.1f);
}