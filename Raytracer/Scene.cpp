#define _USE_MATH_DEFINES
#define MAX_RECURSION_DEPTH 5

#include "Scene.h"

#include <Intersections.h>
#include <math.h>
#include <ObjReader.h>


Scene::Scene() :
  useOctree_(false)
{
  // Add lights
  //lights.push_back(Light(Vector3D(0.5, 5.0, 7.0), ColorD(10.0, 10.0, 10.0)));
  //lights.push_back(Light(Vector3D(0.5, -5.0, 4.0), ColorD(10.0, 0.0, 0.0)));
}

Scene::~Scene()
{

}

bool Scene::Render(uchar* imageData, bool useOctree, int minTriangles, int maxDepth)
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

  double tanHalfFovY = camera.ZNear() * tan(0.5 * camera.FovY());
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

      cameraRays.push_back(cameraRay);

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

  return true;
}

//--------------------------------------------------------------------------------
ColorD Scene::traceRay(Ray ray, double refractiveIndex, int recursionDepth) const
{
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
          hitTime = t;
          hitTriangle = tri;
          hit = true;
        }
      }
    }
  }

  if (hit)
    return radiance(Intersection(ray, hitTime, hitTriangle), ray, refractiveIndex, --recursionDepth);
  else
    return ray.Color * ColorD();  // Background is black
}

//--------------------------------------------------------------------------------
ColorD Scene::radiance(const Intersection& intersection, Ray ray, double refractiveIndex, int recursionDepth) const
{
  ColorD total;

  if (recursionDepth > 0 && ray.Color.IsSignificant())
  {
    double opacity = 1.0;
    double transparency = 1.0 - opacity;

    // Diffuse
    if (opacity > 0.0)
    {
      ColorD surfaceColor = intersection.hit.surfaceColor(intersection.hitPoint);
      total += opacity * surfaceColor * ray.Color * calculateDiffuse(intersection);
    }

    // Reflection
    if (transparency > 0.0)
    {
    }

    // Refraction
  }

  return total;
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
ColorD Scene::calculateReflection() const
{
  ColorD total;

  return total;
}

//--------------------------------------------------------------------------------
ColorD Scene::calculateRefraction() const
{
  ColorD total;

  return total;
}

//--------------------------------------------------------------------------------
void Scene::LoadDefaultScene()
{
	ObjReader reader;
	objects.clear();
  lights.clear();

	Object obj = Object(reader.parseFile("sphere.obj"), Material(ReflectionType::diffuse, ColorD(), ColorD(), 1, 0));
  for (int i = 0; i < obj.triangles.size(); ++i)
  {
    for (int j = 0; j < 3; j++)
    {
      obj.triangles[i].Vertices[j].Position /= 3;
      obj.triangles[i].Vertices[j].Position.X -= 0.5;
      obj.triangles[i].Vertices[j].Position.Z += 0.125;
    }
  }
	objects.push_back(obj);

	obj = Object(reader.parseFile("sphere.obj"), Material(ReflectionType::diffuse, ColorD(), ColorD(), 1, 0));
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
  lights.push_back(Light(Vector3D(-3.0, -5.0, -4.0), ColorD(10.0, 10.0, 10.0)));
  lights.push_back(Light(Vector3D(3.0, 5.0, 4.0), ColorD(10.0, 10.0, 10.0)));

	camera = Camera(Vector3D(0, 0.5, -2.5), Vector3D::Normalise(Vector3D(0, -0.25, 1)), Vector3D(0, 1, 0));
}