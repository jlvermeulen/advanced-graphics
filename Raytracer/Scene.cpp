#define _USE_MATH_DEFINES
#define MAX_RECURSION_DEPTH 5

#include "Scene.h"

#include <Intersections.h>
#include <math.h>
#include <ObjReader.h>
#include <random>


Scene::Scene() :
  useOctree_(false)
{
	// checkerboard
	Vertex v1(Vector3D(-100, -0.5, -100), Vector3D(0, 1, 0), ColorD(), Vector3D());
	Vertex v2(Vector3D(100, -0.5, -100), Vector3D(0, 1, 0), ColorD(), Vector3D());
	Vertex v3(Vector3D(-100, -0.5, 100), Vector3D(0, 1, 0), ColorD(), Vector3D());
	Vertex v4(Vector3D(100, -0.5, 100), Vector3D(0, 1, 0), ColorD(), Vector3D());

	std::deque<Triangle> tris;
	tris.push_back(Triangle(v1, v3, v2));
	tris.push_back(Triangle(v4, v2, v3));

	checkerboard = Object(tris, Material(ReflectionType::diffuse, ColorD(), ColorD(), 1, 0));
}

Scene::~Scene()
{

}

bool Scene::Render(uchar* imageData, bool useOctree, int minTriangles, int maxDepth, RayDistributionType distribution, int numberOfRays, double sigma, int stratificationSize)
{
  useOctree_ = useOctree;

  // Instantiate octrees
  if (useOctree_)
  {
    for (Object& obj : objects)
    {
      obj.CreateOctree(minTriangles, maxDepth);
    }
  }

  switch (distribution)
  {
    case RayDistributionType::none:
      normalRayTrace(imageData);
      break;

    case RayDistributionType::gaussian:
      gaussianRayTrace(imageData, numberOfRays, sigma);
      break;

    case RayDistributionType::jitteredStratification:
      jitteredStratificationRayTrace(imageData, numberOfRays, stratificationSize);
      break;

    case RayDistributionType::stratification:
      stratificationRayTrace(imageData, numberOfRays, stratificationSize);
      break;

    case RayDistributionType::uniform:
      uniformRayTrace(imageData, numberOfRays);
      break;
  }

  return true;
}

void Scene::normalRayTrace(uchar* imageData)
{
  double tanHalfFovY = tan(camera.FovY() / 360 * M_PI);
  double tanHalfFovX = tanHalfFovY * camera.Width / camera.Height;

  double left = -tanHalfFovX;
  double right = tanHalfFovX;
  double top = tanHalfFovY;
  double bottom = -tanHalfFovY;

  // Calculate pixel rays
  for (int x = 0; x < camera.Width; ++x)
  {
    double a = left + (right - left) * (x + 0.5) / camera.Width;
    for (int y = 0; y < camera.Height; ++y)
    {
      double b = top + (bottom - top) * (y + 0.5) / camera.Height;
      Vector3D direction = camera.Focus() + a * camera.Right() + b * camera.Up();

      ColorD intensity(1.0, 1.0, 1.0);
      Ray cameraRay(camera.Eye(), direction, intensity);

      ColorD color = traceRay(cameraRay, 1.0, MAX_RECURSION_DEPTH);

      // clamp color value between 0 and 1
      color.R = std::max(0.0, std::min(1.0, color.R));
      color.G = std::max(0.0, std::min(1.0, color.G));
      color.B = std::max(0.0, std::min(1.0, color.B));

      int offset = (y * camera.Width + x) * 4;

      // For each color channel in reversed order (i.e. blue-green-red-alpha)
      imageData[offset] = (uchar) (color.B * 255.0);
      imageData[offset + 1] = (uchar) (color.G * 255.0);
      imageData[offset + 2] = (uchar) (color.R * 255.0);
      imageData[offset + 3] = 255;
    }
  }
}

//--------------------------------------------------------------------------------
void Scene::gaussianRayTrace(uchar* imageData, int numberOfRays, double sigma)
{
  double tanHalfFovY = tan(camera.FovY() / 360 * M_PI);
  double tanHalfFovX = tanHalfFovY * camera.Width / camera.Height;

  double left = -tanHalfFovX;
  double right = tanHalfFovX;
  double top = tanHalfFovY;
  double bottom = -tanHalfFovY;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dis(0, 1);

  // Calculate pixel rays
  for (int x = 0; x < camera.Width; ++x)
  {
    for (int y = 0; y < camera.Height; ++y)
    {
      ColorD color;
      double totalWeight = 0.0;

      for (int r = 0; r < numberOfRays; ++r)
      {
        double rX = dis(gen);
        double rY = dis(gen);

        double weight = gaussianWeight(rX, rY, sigma);

        double a = left + (right - left) * (x + rX) / camera.Width;
        double b = top + (bottom - top) * (y + rY) / camera.Height;
        Vector3D direction = camera.Focus() + a * camera.Right() + b * camera.Up();

        ColorD intensity(1.0, 1.0, 1.0);
        Ray cameraRay(camera.Eye(), direction, intensity);

        color += weight * traceRay(cameraRay, 1.0, MAX_RECURSION_DEPTH);
        totalWeight += weight;
      }

      color /= totalWeight;

      // clamp color value between 0 and 1
      color.R = std::max(0.0, std::min(1.0, color.R));
      color.G = std::max(0.0, std::min(1.0, color.G));
      color.B = std::max(0.0, std::min(1.0, color.B));

      int offset = (y * camera.Width + x) * 4;

      // For each color channel in reversed order (i.e. blue-green-red-alpha)
      imageData[offset] = (uchar) (color.B * 255.0);
      imageData[offset + 1] = (uchar) (color.G * 255.0);
      imageData[offset + 2] = (uchar) (color.R * 255.0);
      imageData[offset + 3] = 255;
    }
  }
}

//--------------------------------------------------------------------------------
void Scene::jitteredStratificationRayTrace(uchar* imageData, int numberOfRays, int size)
{
}

//--------------------------------------------------------------------------------
void Scene::stratificationRayTrace(uchar* imageData, int numberOfRays, int size)
{
  double tanHalfFovY = tan(camera.FovY() / 360 * M_PI);
  double tanHalfFovX = tanHalfFovY * camera.Width / camera.Height;

  double left = -tanHalfFovX;
  double right = tanHalfFovX;
  double top = tanHalfFovY;
  double bottom = -tanHalfFovY;

  // Calculate pixel rays
  for (int x = 0; x < camera.Width; ++x)
  {
    for (int y = 0; y < camera.Height; ++y)
    {
      ColorD color;

      for (int u = 0; u < size; ++u)
      {
        double mU = (u + 0.5) / size;

        double a = left + (right - left) * (x + mU) / camera.Width;

        for (int v = 0; u < size; ++v)
        {
          double mV = (v + 0.5) / size;

          double b = top + (bottom - top) * (y + mV) / camera.Height;
          Vector3D direction = camera.Focus() + a * camera.Right() + b * camera.Up();

          ColorD intensity(1.0, 1.0, 1.0);
          Ray cameraRay(camera.Eye(), direction, intensity);

          color += traceRay(cameraRay, 1.0, MAX_RECURSION_DEPTH);
        }
      }

      color /= size * size;

      // clamp color value between 0 and 1
      color.R = std::max(0.0, std::min(1.0, color.R));
      color.G = std::max(0.0, std::min(1.0, color.G));
      color.B = std::max(0.0, std::min(1.0, color.B));

      int offset = (y * camera.Width + x) * 4;

      // For each color channel in reversed order (i.e. blue-green-red-alpha)
      imageData[offset] = (uchar) (color.B * 255.0);
      imageData[offset + 1] = (uchar) (color.G * 255.0);
      imageData[offset + 2] = (uchar) (color.R * 255.0);
      imageData[offset + 3] = 255;
    }
  }
}

//--------------------------------------------------------------------------------
void Scene::uniformRayTrace(uchar* imageData, int numberOfRays)
{
  double tanHalfFovY = tan(camera.FovY() / 360 * M_PI);
  double tanHalfFovX = tanHalfFovY * camera.Width / camera.Height;

  double left = -tanHalfFovX;
  double right = tanHalfFovX;
  double top = tanHalfFovY;
  double bottom = -tanHalfFovY;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dis(0, 1);

  // Calculate pixel rays
  for (int x = 0; x < camera.Width; ++x)
  {
    for (int y = 0; y < camera.Height; ++y)
    {
      ColorD color;

      for (int r = 0; r < numberOfRays; ++r)
      {
        double rX = dis(gen);
        double rY = dis(gen);

        double a = left + (right - left) * (x + rX) / camera.Width;
        double b = top + (bottom - top) * (y + rY) / camera.Height;
        Vector3D direction = camera.Focus() + a * camera.Right() + b * camera.Up();

        ColorD intensity(1.0, 1.0, 1.0);
        Ray cameraRay(camera.Eye(), direction, intensity);

        color += traceRay(cameraRay, 1.0, MAX_RECURSION_DEPTH);
      }

      color /= numberOfRays;

      // clamp color value between 0 and 1
      color.R = std::max(0.0, std::min(1.0, color.R));
      color.G = std::max(0.0, std::min(1.0, color.G));
      color.B = std::max(0.0, std::min(1.0, color.B));

      int offset = (y * camera.Width + x) * 4;

      // For each color channel in reversed order (i.e. blue-green-red-alpha)
      imageData[offset] = (uchar) (color.B * 255.0);
      imageData[offset + 1] = (uchar) (color.G * 255.0);
      imageData[offset + 2] = (uchar) (color.R * 255.0);
      imageData[offset + 3] = 255;
    }
  }
}

//--------------------------------------------------------------------------------
double Scene::gaussianWeight(double x, double y, double sigma) const
{
  double deltaX = x - 0.5;
  double deltaY = y - 0.5;

  return exp(-1.0 / 2 * sigma * sigma * (deltaX * deltaX + deltaY * deltaY));
}

//--------------------------------------------------------------------------------
ColorD Scene::traceRay(Ray ray, double refractiveIndex, int recursionDepth) const
{
  Material hitMaterial;
  Triangle hitTriangle;
  double hitTime = std::numeric_limits<double>::max();
  bool hit = false;

  for (const Object& obj : objects)
  {
    if (useOctree_)
    {
      Triangle tri;
      double t = std::numeric_limits<double>::max();

      if (obj.octree.Query(ray, tri, t) && t < hitTime)
      {
        hitMaterial = obj.material;
        hitTriangle = tri;
        hitTime = t;
        hit = true;
      }
    }
    else
    {
      for (const Triangle& tri : obj.triangles)
      {
        double t = std::numeric_limits<double>::max();

        if (Intersects(ray, tri, t) && t < hitTime)
        {
          hitMaterial = obj.material;
          hitTime = t;
          hitTriangle = tri;
          hit = true;
        }
      }
    }
  }

  // checkerboard
  for (const Triangle& tri : checkerboard.triangles)
  {
	double t = std::numeric_limits<double>::max();

	if (Intersects(ray, tri, t) && t < hitTime)
	{
		Vector3D point = ray.Origin + t * ray.Direction;
		double intpart;
		bool xOdd = point.X > 0 && modf(point.X, &intpart) < 0.5 || modf(point.X, &intpart) < -0.5;
		bool zOdd = point.Z > 0 && modf(point.Z, &intpart) < 0.5 || modf(point.Z, &intpart) < -0.5;
		return ray.Color * (xOdd && zOdd || !xOdd && !zOdd ? ColorD(0.0, 0.0, 0.0) : ColorD(1.0, 1.0, 1.0));
	}
  }

  if (hit)
    return radiance(Intersection(ray, hitTime, hitTriangle, hitMaterial), ray, refractiveIndex, --recursionDepth);
  else
    return ray.Color * ColorD();  // Background is black
}

//--------------------------------------------------------------------------------
ColorD Scene::radiance(const Intersection& intersection, Ray ray, double refractiveIndex, int recursionDepth) const
{
  if (recursionDepth > 0 && ray.Color.IsSignificant())
  {
    double transparency = intersection.hitMaterial.transparency;
    double opacity = 1.0 - transparency;

    // Diffuse
    if (intersection.hitMaterial.reflType == ReflectionType::diffuse)
    {
      if (opacity > 0.0)
      {
        ColorD surfaceColor = intersection.hit.surfaceColor(intersection.hitPoint);
        return opacity * surfaceColor * ray.Color * calculateDiffuse(intersection);
      }
    }

    // Refraction
    if (intersection.hitMaterial.reflType == ReflectionType::refractive)
    {
      if (transparency > 0.0)
      {
        return transparency * calculateRefraction(intersection, ray, refractiveIndex, recursionDepth);
      }
    }

    // Reflection
    if (intersection.hitMaterial.reflType == ReflectionType::specular)
    {
      return opacity * calculateReflection(intersection, ray, refractiveIndex, recursionDepth);
    }
  }

  return ColorD();
}

//--------------------------------------------------------------------------------
ColorD Scene::calculateDiffuse(const Intersection& intersection) const
{
  ColorD total;

  for (const Light& source : lights)
  {
    Ray toSource(intersection.hitPoint, source.position - intersection.hitPoint, ColorD());

    Triangle blockingTriangle;
    double blockingTime;

    // Check for occluding object triangles
    // between the hitPoint and the light source
    bool noOcclusion = true;

    for (const Object& obj : objects)
    {
      if (useOctree_)
      {
        if (obj.octree.Query(toSource, blockingTriangle, blockingTime) && blockingTime < 1.0)
        {
          noOcclusion = false;
          break;
        }
      }
      else
      {
        for (const Triangle& tri : obj.triangles)
        {
          if (Intersects(toSource, tri, blockingTime) && blockingTime < 1.0)
          {
            noOcclusion = false;
            break;
          }
        }
      }
    }

    // No triangles occluding the light source
    if (noOcclusion)
    {
      double incidence = Vector3D::Dot(intersection.hit.surfaceNormal(intersection.hitPoint), Vector3D::Normalise(toSource.Direction));

      if (incidence > 0.0)
      {
        double intensity = incidence / toSource.Direction.LengthSquared();

        total += intensity * source.color;
      }
    }
  }

  return total;
}

//--------------------------------------------------------------------------------
ColorD Scene::calculateReflection(const Intersection& intersection, const Ray& ray, double refractiveIndex, int recursionDepth) const
{
  const Vector3D& normal = intersection.hit.surfaceNormal(intersection.hitPoint);
  Vector3D reflectDir = ray.Direction - (2 * Vector3D::Dot(ray.Direction, normal)) * normal;

  Ray out(intersection.hitPoint, reflectDir, ray.Color);

  return traceRay(out, refractiveIndex, --recursionDepth);
}

//--------------------------------------------------------------------------------
ColorD Scene::calculateRefraction(const Intersection& intersection, const Ray& ray, double refractiveIndex, int recursionDepth) const
{
  Vector3D normal = intersection.hit.surfaceNormal(intersection.hitPoint);

  double c = Vector3D::Dot(-normal, ray.Direction);
  double r = refractiveIndex / intersection.hitMaterial.refrIndex;
  double rad = 1.0 - r * r * (1.0 - c * c);

  if (rad < 0)
    return calculateReflection(intersection, ray, refractiveIndex, recursionDepth);

  Vector3D refractDir = r * ray.Direction + (r * c - sqrt(rad)) * normal;

  Ray out(intersection.hitPoint, refractDir, ray.Color);

  return traceRay(out, intersection.hitMaterial.refrIndex, --recursionDepth);
}

//--------------------------------------------------------------------------------
void Scene::LoadDefaultScene()
{
	ObjReader reader;
	objects.clear();
  lights.clear();

	Object obj = Object(reader.parseFile("sphere.obj"), Material(ReflectionType::diffuse, ColorD(), ColorD(), 1.0, 0.0));
  for (int i = 0; i < obj.triangles.size(); ++i)
  {
    for (int j = 0; j < 3; j++)
    {
      obj.triangles[i].Vertices[j].Position /= 3;
      obj.triangles[i].Vertices[j].Position.X -= 0.5;
      obj.triangles[i].Vertices[j].Position.Z += 0.125;
      obj.triangles[i].Vertices[j].Color = ColorD(1.0, 0.0, 0.0);
    }
  }
	objects.push_back(obj);

	obj = Object(reader.parseFile("sphere.obj"), Material(ReflectionType::specular, ColorD(), ColorD(), 0.5, 0.5));
  for (int i = 0; i < obj.triangles.size(); ++i)
  {
    for (int j = 0; j < 3; j++)
    {
      obj.triangles[i].Vertices[j].Position /= 3;
      obj.triangles[i].Vertices[j].Position.X += 0.5;
      obj.triangles[i].Vertices[j].Position.Z -= 0.125;
    }
  }
	objects.push_back(obj);

  // Add lights
  lights.push_back(Light(Vector3D(-3.0, -5.0, -4.0), ColorD(15.0, 15.0, 15.0)));
  lights.push_back(Light(Vector3D(3.0, 5.0, 4.0), ColorD(15.0, 15.0, 15.0)));

	camera = Camera(Vector3D(0, 15.0, -10.0), Vector3D::Normalise(Vector3D(0, -0.5, 1)), Vector3D(0, 1, 0));
}