#define _USE_MATH_DEFINES
#define MAX_RECURSION_DEPTH 5

#include "Scene.h"

#include <Intersections.h>
#include <math.h>
#include <ObjReader.h>

Scene::Scene() :
  useOctree_(false)
{
	// checkerboard
	Vertex v1(Vector3D(-100, -0.5, -100), Vector3D(0, 1, 0), ColorD(1.0, 1.0, 1.0), Vector3D());
	Vertex v2(Vector3D(100, -0.5, -100), Vector3D(0, 1, 0), ColorD(1.0, 1.0, 1.0), Vector3D());
	Vertex v3(Vector3D(-100, -0.5, 100), Vector3D(0, 1, 0), ColorD(1.0, 1.0, 1.0), Vector3D());
	Vertex v4(Vector3D(100, -0.5, 100), Vector3D(0, 1, 0), ColorD(1.0, 1.0, 1.0), Vector3D());

	std::deque<Triangle> tris;
	tris.push_back(Triangle(v1, v3, v2));
	tris.push_back(Triangle(v4, v2, v3));

	//checkerboard = Object(tris, Material(ReflectionType::diffuse, ColorD(1.0, 1.0, 1.0), ColorD(), 1, 0));

	// Add lights
	lights.push_back(Light(Vector3D(-3.0, -5.0, -4.0), ColorD(25.0, 25.0, 25.0)));
	lights.push_back(Light(Vector3D(3.0, 5.0, 4.0), ColorD(25.0, 25.0, 25.0)));

	camera = Camera(Vector3D(0, 15.0, -10.0), Vector3D::Normalise(Vector3D(0, -0.5, 1)), Vector3D(0, 1, 0));
}

Scene::~Scene()
{

}

bool Scene::Render(uchar* imageData, bool useOctree, int minTriangles, int maxDepth, int samplesPerPixel, double sigma)
{
	// Instantiate octrees
	for (Object& obj : objects)
		obj.CreateOctree(minTriangles, maxDepth);

	std::pair<ColorD, ColorD>* samples = new std::pair<ColorD, ColorD>[camera.Width * camera.Height];
	tracePixels(samples, samplesPerPixel, sigma);

	#pragma omp parallel for
	for (int x = 0; x < camera.Width; ++x)
		for (int y = 0; y < camera.Height; ++y)
		{
			std::pair<ColorD, ColorD> sample = samples[y * camera.Width + x];
			ColorD color = sample.first / sample.second;
			color.Clip();

			int offset = (y * camera.Width + x) * 4;

			// For each color channel in reversed order (i.e. blue-green-red-alpha)
			imageData[offset] = (uchar) (color.B * 255.0);
			imageData[offset + 1] = (uchar) (color.G * 255.0);
			imageData[offset + 2] = (uchar) (color.R * 255.0);
			imageData[offset + 3] = 255;
		}

	delete [] samples;

	return true;
}

//--------------------------------------------------------------------------------
void Scene::tracePixels(std::pair<ColorD, ColorD>* pixelData, int samplesPerPixel, double sigma)
{
	double tanHalfFovY = tan(camera.FovY() / 360 * M_PI);
	double tanHalfFovX = tanHalfFovY * camera.Width / camera.Height;

	double left = -tanHalfFovX;
	double right = tanHalfFovX;
	double top = tanHalfFovY;
	double bottom = -tanHalfFovY;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dist(0, 1);

	// Calculate pixel rays
	for (int a = 0; a < 3; a++) // prevent concurrent access of same array indices
	{
		#pragma omp parallel for schedule(dynamic)
		for (int x = a; x < camera.Width; x += 3)
			for (int y = 0; y < camera.Height; ++y)
			{
				for (int r = 0; r < samplesPerPixel; ++r)
				{
					double rX = dist(gen);
					double rY = dist(gen);

					double a = left + (right - left) * (x + rX) / camera.Width;
					double b = top + (bottom - top) * (y + rY) / camera.Height;

					Vector3D direction = Vector3D::Normalise(camera.Focus() + a * camera.Right() + b * camera.Up());
					Ray cameraRay(camera.Eye(), direction);

					for (unsigned int c = 0; c < 3; ++c)
					{
						double value = traceRay(cameraRay, c, dist, gen);

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

								double weight = gaussianWeight(i - rX - 0.5, j - rY - 0.5, sigma);

								int index = yy * camera.Width + xx;
								std::pair<ColorD, ColorD>& data = pixelData[index];

								data.first[c] += value * weight;
								data.second[c] += weight;
							}
						}
					}
				}
			}
	}
}

//--------------------------------------------------------------------------------
double Scene::gaussianWeight(double dx, double dy, double sigma) const
{
  return exp(-(dx * dx + dy * dy) / (2 * sigma * sigma));
}

//--------------------------------------------------------------------------------
double Scene::traceRay(const Ray& ray, unsigned int channel, const std::uniform_real_distribution<double>& dist, std::mt19937& gen) const
{
	Material hitMaterial;
	Triangle hitTriangle;
	double hitTime = std::numeric_limits<double>::max();
	bool hit = false;
	for (const Object& obj : objects)
	{
		Triangle tri;
		double t = std::numeric_limits<double>::max();

		if (obj.octree->Query(ray, tri, t) && t < hitTime)
		{
			hitMaterial = obj.material;
			hitTriangle = tri;
			hitTime = t;
			hit = true;
		}
	}

	// checkerboard
	/*for (const Triangle& tri : checkerboard.triangles)
	{
		double t = std::numeric_limits<double>::max();

		if (Intersects(ray, tri, t) && t < hitTime)
		{
			Vector3D point = ray.Origin + t * ray.Direction;
			double intpart;
			bool xOdd = point.X > 0 && modf(point.X, &intpart) < 0.5 || modf(point.X, &intpart) < -0.5;
			bool zOdd = point.Z > 0 && modf(point.Z, &intpart) < 0.5 || modf(point.Z, &intpart) < -0.5;
			hitTime = t;
			hitTriangle = tri;
			hit = xOdd && zOdd || !xOdd && !zOdd;
		}
	}*/

	if (hit)
	{
		if (hitMaterial.emission.IsSignificant()) // hit a light source
			return hitMaterial.emission[channel];

		if (dist(gen) > hitMaterial.color[channel]) // russian roulette
			return 0;

		const Vector3D& hitPoint = ray.Origin + hitTime * ray.Direction;
		const Vector3D& normal = hitTriangle.surfaceNormal(hitPoint);

		if (ray.Direction.Dot(normal) > 0) // hit from behind
			return 0;

		if (hitMaterial.reflType == ReflectionType::specular)
			return traceRay(Ray(hitPoint, Vector3D::Normalise(Vector3D::Reflect(ray.Direction, normal))), channel, dist, gen); // perfect reflection
		else if (hitMaterial.reflType == ReflectionType::diffuse)
		{
			double u = dist(gen) * 2 - 1, theta = dist(gen) * M_PI * 2, x = sqrt(1 - u * u); // sample unit sphere
			Vector3D hemi(x * cos(theta), x * sin(theta), u);
			if (hemi.Dot(normal) < 0) // flip if "behind" normal
				hemi *= -1;

			return traceRay(Ray(hitPoint, hemi), channel, dist, gen);
		}
		else if (hitMaterial.reflType == ReflectionType::refractive)
			return traceRay(Ray(hitPoint, Vector3D::Normalise(Vector3D::Refract(ray.Direction, normal, 1.0, hitMaterial.refrIndex))), channel, dist, gen);
		else
			return 0;
	}
	else
		return 0;  // Background is black
}

//--------------------------------------------------------------------------------
//ColorD Scene::radiance(const Intersection& intersection, Ray ray, double refractiveIndex, int recursionDepth) const
//{
//  if (recursionDepth > 0 && ray.Color.IsSignificant())
//  {
//    double transparency = intersection.hitMaterial.transparency;
//    double opacity = 1.0 - transparency;
//
//    // Diffuse
//    if (intersection.hitMaterial.reflType == ReflectionType::diffuse)
//    {
//      if (opacity > 0.0)
//      {
//        ColorD surfaceColor = intersection.hit.surfaceColor(intersection.hitPoint);
//        return opacity * surfaceColor * ray.Color * calculateDiffuse(intersection);
//      }
//    }
//
//    // Refraction
//    if (intersection.hitMaterial.reflType == ReflectionType::refractive)
//    {
//      if (transparency > 0.0)
//      {
//        return transparency * calculateRefraction(intersection, ray, refractiveIndex, recursionDepth);
//      }
//    }
//
//    // Reflection
//    if (intersection.hitMaterial.reflType == ReflectionType::specular)
//    {
//      return opacity * calculateReflection(intersection, ray, refractiveIndex, recursionDepth);
//    }
//  }
//
//  return ColorD();
//}

//--------------------------------------------------------------------------------
//ColorD Scene::calculateDiffuse(const Intersection& intersection) const
//{
//  ColorD total;
//
//  for (const Light& source : lights)
//  {
//    Ray toSource(intersection.hitPoint, source.position - intersection.hitPoint, ColorD());
//
//    Triangle blockingTriangle;
//    double blockingTime;
//
//    // Check for occluding object triangles
//    // between the hitPoint and the light source
//    bool noOcclusion = true;
//
//    for (const Object& obj : objects)
//    {
//      if (useOctree_)
//      {
//        if (obj.octree->Query(toSource, blockingTriangle, blockingTime) && blockingTime < 1.0)
//        {
//          noOcclusion = false;
//          break;
//        }
//      }
//      else
//      {
//        for (const Triangle& tri : obj.triangles)
//        {
//          if (Intersects(toSource, tri, blockingTime) && blockingTime < 1.0)
//          {
//            noOcclusion = false;
//            break;
//          }
//        }
//      }
//    }
//
//    // No triangles occluding the light source
//    if (noOcclusion)
//    {
//      double incidence = Vector3D::Dot(intersection.hit.surfaceNormal(intersection.hitPoint), Vector3D::Normalise(toSource.Direction));
//
//      if (incidence > 0.0)
//      {
//        double intensity = incidence / toSource.Direction.LengthSquared();
//
//        total += intensity * source.color;
//      }
//    }
//  }
//
//  return total;
//}

//--------------------------------------------------------------------------------
//ColorD Scene::calculateReflection(const Intersection& intersection, const Ray& ray, double refractiveIndex, int recursionDepth) const
//{
//  const Vector3D& normal = intersection.hit.surfaceNormal(intersection.hitPoint);
//  Vector3D reflectDir = ray.Direction - (2 * Vector3D::Dot(ray.Direction, normal)) * normal;
//
//  Ray out(intersection.hitPoint, reflectDir, ray.Color);
//
//  return traceRay(out, refractiveIndex, --recursionDepth);
//}

//--------------------------------------------------------------------------------
//ColorD Scene::calculateRefraction(const Intersection& intersection, const Ray& ray, double refractiveIndex, int recursionDepth) const
//{
//  Vector3D normal = intersection.hit.surfaceNormal(intersection.hitPoint);
//
//  double c = Vector3D::Dot(-normal, ray.Direction);
//  double r = refractiveIndex / intersection.hitMaterial.refrIndex;
//  double rad = 1.0 - r * r * (1.0 - c * c);
//
//  if (rad < 0)
//    return calculateReflection(intersection, ray, refractiveIndex, recursionDepth);
//
//  Vector3D refractDir = r * ray.Direction + (r * c - sqrt(rad)) * normal;
//
//  Ray out(intersection.hitPoint, refractDir, ray.Color);
//
//  return traceRay(out, intersection.hitMaterial.refrIndex, --recursionDepth);
//}

//--------------------------------------------------------------------------------
Ray Scene::SampleLight(std::vector<Lightarea> lightsources, Vector3D hitPoint) const
{
	Triangle triangle;
	float emission;
	const int n = lightsources.size();
	std::vector<float> flux;
	float totalflux = 0;
	// loop over all lightsources
	for (int l = 0; l < lightsources.size(); l++)
	{
		triangle = lightsources[l].triangle;
		emission = lightsources[l].emission;
		Vector3D outgoingRay = hitPoint - triangle.center();
		float distance = outgoingRay.Length();
		Vector3D normal = triangle.surfaceNormal(triangle.center());
		normal.Normalise();
		outgoingRay.Normalise();
		
		flux.push_back(triangle.area() * emission * std::abs(Vector3D::Dot(outgoingRay, normal)) / distance);
		totalflux += flux[l];
	}
	// devide all flux values by the total flux gives the probabilities
	for (int l = 0; l < lightsources.size(); l++)
		flux[l] /= totalflux;

	// get the light using the probabilities and a random float between 0 and 1
	float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
	int light = -1;
	for (int l = 0; l < lightsources.size(); l++)
	{
		if (flux[l] >= r)
		{
			light = l;
			break;
		}
		else
			r -= flux[l];
	}
	// ensure that we have a light selected (might not be because of float roundings)
	if (light == -1) light = lightsources.size() - 1;

	// create random point withing the lightsource and return ray towards it.
	float u = 1;
	float v = 1;
	while (u + v > 1)
	{
		u = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		v = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
	}
	triangle = lightsources[light].triangle;
	Vector3D v1 = triangle.Vertices[1].Position - triangle.Vertices[0].Position;
	Vector3D v2 = triangle.Vertices[2].Position - triangle.Vertices[0].Position;
	return Ray(hitPoint, (v1*u + v2*v) - hitPoint);


}

//--------------------------------------------------------------------------------
void Scene::LoadDefaultScene()
{
  ObjReader reader;
  objects.clear();
  lights.clear();

  // left sphere
  Object obj = Object(reader.parseFile("sphere.obj"), Material(ReflectionType::diffuse, ColorD(1.0, 0.5, 0.0), ColorD(), 1.5, 0.0));
  for (unsigned int i = 0; i < obj.triangles.size(); ++i)
  {
    for (int j = 0; j < 3; j++)
    {
      obj.triangles[i].Vertices[j].Position /= 3;
      obj.triangles[i].Vertices[j].Position.X -= 0.5;
      obj.triangles[i].Vertices[j].Position.Z += 0.125;
      obj.triangles[i].Vertices[j].Color = obj.material.color;
    }
  }
  objects.push_back(obj);

  // right sphere
  obj = Object(reader.parseFile("sphere.obj"), Material(ReflectionType::specular, ColorD(1.0, 1.0, 1.0), ColorD(), 0.5, 0.5));
  for (unsigned int i = 0; i < obj.triangles.size(); ++i)
  {
    for (int j = 0; j < 3; j++)
    {
      obj.triangles[i].Vertices[j].Position /= 3;
      obj.triangles[i].Vertices[j].Position.X += 0.5;
      obj.triangles[i].Vertices[j].Position.Z -= 0.125;
	  obj.triangles[i].Vertices[j].Color = obj.material.color;
    }
  }
  objects.push_back(obj);

  // light
  obj = Object(reader.parseFile("sphere.obj"), Material(ReflectionType::specular, ColorD(1.0, 1.0, 1.0), ColorD(25.0, 25.0, 25.0), 0.5, 0.5));
  for (unsigned int i = 0; i < obj.triangles.size(); ++i)
  {
    for (int j = 0; j < 3; j++)
    {
      obj.triangles[i].Vertices[j].Position /= 1.5;
	  obj.triangles[i].Vertices[j].Position.Y /= 12;
      obj.triangles[i].Vertices[j].Position.Y += 1.5;
	  obj.triangles[i].Vertices[j].Color = obj.material.color;
    }
  }
  objects.push_back(obj);

  // right
  obj = Object(reader.parseFile("cube.obj"), Material(ReflectionType::diffuse, ColorD(0.5, 0, 0), ColorD(), 1.0, 0.0));
  for (unsigned int i = 0; i < obj.triangles.size(); ++i)
  {
    for (int j = 0; j < 3; j++)
    {
      obj.triangles[i].Vertices[j].Position *= 5;
      obj.triangles[i].Vertices[j].Position.X += 4;
	  obj.triangles[i].Vertices[j].Color = obj.material.color;
    }
  }
  objects.push_back(obj);

  // left
  obj = Object(reader.parseFile("cube.obj"), Material(ReflectionType::diffuse, ColorD(0, 0, 0.5), ColorD(), 1.0, 0.0));
  for (unsigned int i = 0; i < obj.triangles.size(); ++i)
  {
    for (int j = 0; j < 3; j++)
    {
      obj.triangles[i].Vertices[j].Position *= 5;
      obj.triangles[i].Vertices[j].Position.X -= 4;
	  obj.triangles[i].Vertices[j].Color = obj.material.color;
    }
  }
  objects.push_back(obj);

  // back
  obj = Object(reader.parseFile("cube.obj"), Material(ReflectionType::diffuse, ColorD(0.5, 0.5, 0.5), ColorD(), 1.0, 0.0));
  for (unsigned int i = 0; i < obj.triangles.size(); ++i)
  {
    for (int j = 0; j < 3; j++)
    {
      obj.triangles[i].Vertices[j].Position *= 5;
      obj.triangles[i].Vertices[j].Position.Z -= 4;
	  obj.triangles[i].Vertices[j].Color = obj.material.color;
    }
  }
  objects.push_back(obj);

  // top
  obj = Object(reader.parseFile("cube.obj"), Material(ReflectionType::diffuse, ColorD(0.0, 0.5, 0.0), ColorD(), 1.0, 0.0));
  for (unsigned int i = 0; i < obj.triangles.size(); ++i)
  {
    for (int j = 0; j < 3; j++)
    {
      obj.triangles[i].Vertices[j].Position *= 5;
      obj.triangles[i].Vertices[j].Position.Y += 4;
	  obj.triangles[i].Vertices[j].Color = obj.material.color;
    }
  }
  objects.push_back(obj);

  // bottom
  obj = Object(reader.parseFile("cube.obj"), Material(ReflectionType::diffuse, ColorD(0.5, 0.5, 0.0), ColorD(), 1.0, 0.0));
  for (unsigned int i = 0; i < obj.triangles.size(); ++i)
  {
    for (int j = 0; j < 3; j++)
    {
      obj.triangles[i].Vertices[j].Position *= 5;
      obj.triangles[i].Vertices[j].Position.Y -= 4;
	  obj.triangles[i].Vertices[j].Color = obj.material.color;
    }
  }
  objects.push_back(obj);

  // Add lights
  lights.push_back(Light(Vector3D(-0.75, 1.25, -1.0), ColorD(5.0, 5.0, 5.0)));
  lights.push_back(Light(Vector3D(0.75, 1.25, 1.0), ColorD(5.0, 5.0, 5.0)));

  // Add lightAreas
 /* Vertex v1 = Vertex(Vector3D(-0.6, 1.25, 1.0), Vector3D(0, 1, 0), ColorD(1.0, 1.0, 1.0), Vector3D(0, 0, 0));
  Vertex v2 = Vertex(Vector3D(-0.8, 1.25, 1.0), Vector3D(0, 1, 0), ColorD(1.0, 1.0, 1.0), Vector3D(0, 0, 0));
  Vertex v3 = Vertex(Vector3D(-0.7, 1.25, 1.1), Vector3D(0, 1, 0), ColorD(1.0, 1.0, 1.0), Vector3D(0, 0, 0));
  Triangle t1 = Triangle(v1, v2, v3);
  Vertex v12 = Vertex(Vector3D(0.6, 1.25, 1.0), Vector3D(0, 1, 0), ColorD(1.0, 1.0, 1.0), Vector3D(0, 0, 0));
  Vertex v22 = Vertex(Vector3D(0.8, 1.25, 1.0), Vector3D(0, 1, 0), ColorD(1.0, 1.0, 1.0), Vector3D(0, 0, 0));
  Vertex v32 = Vertex(Vector3D(0.7, 1.25, 1.1), Vector3D(0, 1, 0), ColorD(1.0, 1.0, 1.0), Vector3D(0, 0, 0));
  Triangle t2 = Triangle(v12, v22, v32);
  lightareas.push_back(Lightarea(t1, ColorD(1.0, 1.0, 1.0), 10));
  lightareas.push_back(Lightarea(t2, ColorD(1.0, 1.0, 1.0), 10));*/

  //camera = Camera(Vector3D(0, 15.0, -10.0), Vector3D::Normalise(Vector3D(0, -0.5, 1)), Vector3D(0, 1, 0));
  camera = Camera(Vector3D(0, 0, 4.5), Vector3D::Normalise(Vector3D(0, 0, -1)), Vector3D(0, 1, 0));
}