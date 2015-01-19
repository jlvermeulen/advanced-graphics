#define _USE_MATH_DEFINES
#define RUSSIAN_ROULETTE_PROBABILITY 0.90
#define MIN_PATH_LENGTH 2

#include "Scene.h"

#include "Intersections.h"
#include "ObjReader.h"

#include <math.h>
#include <random>

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

bool Scene::Render(uchar* imageData, int minTriangles, int maxDepth, int samplesPerPixel, double sigma)
{
	// Instantiate octrees
	for (Object& obj : objects)
      obj.ConstructOctree(minTriangles, maxDepth);

	std::pair<ColorD, double>* samples = new std::pair<ColorD, double>[camera.Width * camera.Height];
	TracePixels(samples, samplesPerPixel, sigma);

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
void Scene::TracePixels(std::pair<ColorD, double>* pixelData, int samplesPerPixel, double sigma)
{
	double tanHalfFovY = tan(camera.FovY() / 360 * M_PI);
	double tanHalfFovX = tanHalfFovY * camera.Width / camera.Height;

	double left = -tanHalfFovX;
	double right = tanHalfFovX;
	double top = tanHalfFovY;
	double bottom = -tanHalfFovY;

	// Calculate pixel rays
	for (int o = 0; o < 3; ++o) // prevent concurrent access of same array indices
	{
		#pragma omp parallel for schedule(dynamic)
    for (int x = o; x < camera.Width; x += 3)
    {
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
	return material.emission + DirectIllumination(point, normal, material) + IndirectIllumination(point, in, normal, material, depth);
}

ColorD Scene::DirectIllumination(const Vector3D& point, const Vector3D& normal, const Material& material)
{
	std::pair<Ray, Triangle> sample = SampleLight(point);

	Material hitMaterial;
	Triangle hitTriangle;
	double hitTime;

	Ray ray(sample.first.Origin + 0.005 * normal, sample.first.Direction);
	if (!FirstHitInfo(ray, hitTime, hitTriangle, hitMaterial) || hitTriangle != sample.second)
		return ColorD();

	return hitMaterial.emission * material.color;
}

ColorD Scene::IndirectIllumination(Vector3D point, const Vector3D& in, const Vector3D& normal, const Material& material, unsigned int depth)
{
	if (depth > MIN_PATH_LENGTH && dist(gen) > RUSSIAN_ROULETTE_PROBABILITY)
		return ColorD();

	Ray ray(point, in);

	double c = 1 / (2 * M_PI);
  ColorD value(c, c, c);

  // Sample unit sphere
  double phi = 2.0 * M_PI * dist(gen);
  double cosPhi = cos(phi);
  double sinPhi = sin(phi);

  double cosTheta = pow(1.0 - dist(gen), 1.0 / (1.0 + material.specularExponent));
  double sinTheta = sqrt(1.0 - cosTheta * cosTheta);

  double x = sinTheta * cosPhi;
  double y = sinTheta * sinPhi;
  double z = cosTheta;

  Vector3D hemi(x, y, z);

  //double z = dist(gen) * 2 - 1, theta = dist(gen) * M_PI * 2, r = sqrt(1 - z * z);

  //Vector3D hemi(r * cos(theta), r * sin(theta), z);
  //if (hemi.Dot(normal) < 0) // flip if "behind" normal
  //  hemi *= -1;

  if (material.reflType == ReflectionType::specular)
  {
    ray.Reflect(point, normal);

    Vector3D w = ray.Direction;

    Vector3D u = Vector3D::Cross(Vector3D(0.00424, 1, 0.00764), w);
    u.Normalise();

    Vector3D v = Vector3D::Cross(u, w);

    Vector3D reflDir = hemi.X * u + hemi.Y * v + hemi.Z * w;

    if (reflDir.Dot(normal) < 0)
      reflDir *= -1;

    ray = Ray(point, reflDir);
    value *= material.color * pow(w.Dot(reflDir), material.specularExponent) * Vector3D::Dot(normal, reflDir);
  }
	else if (material.reflType == ReflectionType::diffuse)
	{
		ray = Ray(point, hemi);
		point += 0.005 * normal;
		value *= material.color * Vector3D::Dot(normal, ray.Direction) / (2 * M_PI);
	}
	else if (material.reflType == ReflectionType::refractive)
		ray.Refract(point, normal, 1.0, material.refrIndex);

	Material hitMaterial;
	Triangle hitTriangle;
	double hitTime;

	if (!FirstHitInfo(ray, hitTime, hitTriangle, hitMaterial))
		return ColorD();

	Vector3D hitPoint = ray.Origin + hitTime * ray.Direction;
	return ComputeRadiance(hitPoint, ray.Direction, hitTriangle, hitMaterial, depth + 1) * value / RUSSIAN_ROULETTE_PROBABILITY;
}

//--------------------------------------------------------------------------------
bool Scene::FirstHitInfo(const Ray& ray, double& time, Triangle& triangle, Material& mat) const
{
	time = std::numeric_limits<double>::max();
	bool hit = false;

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
			double distance = outgoingRay.Length();
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
  //lights.clear();

  // Example: Two spheres

  // Parse left sphere
  reader.parseFile("../models/sphere.obj");

  // Parse right sphere
  reader.parseFile("../models/sphere.obj");

  // Parse light
  reader.parseFile("../models/light2.obj");

  // Parse light
  reader.parseFile("../models/light3.obj");

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

  // Left sphere
  objects[0].material = Material(ReflectionType::diffuse, ColorD(1.0, 0.5, 0.0), ColorD(), 0.5, 0.0, 0.5);
  unsigned int nTriangles = objects[0].triangles.size();

  for (unsigned int i = 0; i < nTriangles; ++i)
  {
    for (unsigned int j = 0; j < 3; ++j)
    {
      objects[0].triangles[i].Vertices[j].Position /= 3;
      objects[0].triangles[i].Vertices[j].Position.X -= 0.5;
      objects[0].triangles[i].Vertices[j].Position.Z += 0.125;
    }
  }

  // Right sphere
  objects[1].material = Material(ReflectionType::specular, ColorD(0.5, 0.5, 0.5), ColorD(), 0.5, 100.0, 0.5);
  nTriangles = objects[1].triangles.size();

  for (unsigned int i = 0; i < nTriangles; ++i)
  {
    for (unsigned int j = 0; j < 3; ++j)
    {
      objects[1].triangles[i].Vertices[j].Position /= 3;
      objects[1].triangles[i].Vertices[j].Position.X += 0.5;
      objects[1].triangles[i].Vertices[j].Position.Z -= 0.125;
    }
  }

  // Top Light
  objects[2].material = Material(ReflectionType::specular, ColorD(0.0, 0.0, 0.0), ColorD(0.75, 0.75, 0.75), 0.5, 0.0, 0.5);
  nTriangles = objects[2].triangles.size();

  for (unsigned int i = 0; i < nTriangles; ++i)
  {
    for (unsigned int j = 0; j < 3; ++j)
      objects[2].triangles[i].Vertices[j].Position.Y += 1.5;

    objects[2].triangles[i].CalculateArea();
    objects[2].triangles[i].CalculateCenter();
  }

  lights.push_back(objects[2]);

  // Bottom Light
  objects[3].material = Material(ReflectionType::specular, ColorD(0.0, 0.0, 0.0), ColorD(0.4, 0.4, 0.4), 0.5, 0.0, 0.5);
  nTriangles = objects[3].triangles.size();

  for (unsigned int i = 0; i < nTriangles; ++i)
  {
    for (unsigned int j = 0; j < 3; ++j)
      objects[3].triangles[i].Vertices[j].Position.Y -= 1.5;
    objects[3].triangles[i].CalculateArea();
    objects[3].triangles[i].CalculateCenter();
  }

  lights.push_back(objects[3]);

  // Right
  objects[4].material = Material(ReflectionType::diffuse, ColorD(0.9, 0, 0), ColorD(), 1.0, 0.0, 0.0);
  nTriangles = objects[4].triangles.size();

  for (unsigned int i = 0; i < nTriangles; ++i)
  {
    for (unsigned int j = 0; j < 3; ++j)
    {
      objects[4].triangles[i].Vertices[j].Position *= 5;
      objects[4].triangles[i].Vertices[j].Position.X += 4;
    }
  }

  // Left
  objects[5].material = Material(ReflectionType::diffuse, ColorD(0, 0, 0.9), ColorD(), 1.0, 0.0, 0.0);
  nTriangles = objects[5].triangles.size();

  for (unsigned int i = 0; i < nTriangles; ++i)
  {
    for (unsigned int j = 0; j < 3; ++j)
    {
      objects[5].triangles[i].Vertices[j].Position *= 5;
      objects[5].triangles[i].Vertices[j].Position.X -= 4;
    }
  }

  // Back
  objects[6].material = Material(ReflectionType::diffuse, ColorD(0.5, 0.5, 0.5), ColorD(), 1.0, 0.0, 0.0);
  nTriangles = objects[6].triangles.size();

  for (unsigned int i = 0; i < nTriangles; ++i)
  {
    for (unsigned int j = 0; j < 3; ++j)
    {
      objects[6].triangles[i].Vertices[j].Position *= 5;
      objects[6].triangles[i].Vertices[j].Position.Z -= 4;
    }
  }


  // Top
  objects[7].material = Material(ReflectionType::diffuse, ColorD(0.0, 0.9, 0.0), ColorD(), 1.0, 0.0, 0.0);
  nTriangles = objects[7].triangles.size();

  for (unsigned int i = 0; i < nTriangles; ++i)
  {
    for (unsigned int j = 0; j < 3; ++j)
    {
      objects[7].triangles[i].Vertices[j].Position *= 5;
      objects[7].triangles[i].Vertices[j].Position.Y += 4;
    }
  }

  // Bottom
  objects[8].material = Material(ReflectionType::diffuse, ColorD(0.5, 0.5, 0.0), ColorD(), 1.0, 0.0, 0.0);
  nTriangles = objects[8].triangles.size();

  for (unsigned int i = 0; i < nTriangles; ++i)
  {
    for (unsigned int j = 0; j < 3; ++j)
    {
      objects[8].triangles[i].Vertices[j].Position *= 5;
      objects[8].triangles[i].Vertices[j].Position.Y -= 4;
    }
  }

  // Add lights
  //lights.push_back(Light(Vector3D(-0.75, 1.25, -1.0), ColorD(5.0, 5.0, 5.0)));
  //lights.push_back(Light(Vector3D(0.75, 1.25, 1.0), ColorD(5.0, 5.0, 5.0)));

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