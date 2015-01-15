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
	std::uniform_real_distribution<double> dis(0, 1);

	// Calculate pixel rays
	for (int a = 0; a < 3; a++) // prevent concurrent access of same array indices
	{
		#pragma omp parallel for
		for (int x = a; x < camera.Width; x += 3)
			for (int y = 0; y < camera.Height; ++y)
			{
				for (int r = 0; r < samplesPerPixel; ++r)
				{
					double rX = dis(gen);
					double rY = dis(gen);

					double a = left + (right - left) * (x + rX) / camera.Width;
					double b = top + (bottom - top) * (y + rY) / camera.Height;

					Vector3D direction = camera.Focus() + a * camera.Right() + b * camera.Up();
					Ray cameraRay(camera.Eye(), direction);

					for (unsigned int c = 0; c < 3; ++c)
					{
						double value = traceRay(cameraRay, c, MAX_RECURSION_DEPTH);

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
double Scene::traceRay(Ray ray, unsigned int channel, unsigned int recursionDepth) const
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
		// TODO: russian roulette, reflection
		return hitMaterial.color[channel];
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
void Scene::LoadDefaultScene()
{
  ObjReader reader;
  objects.clear();
  lights.clear();

  Object obj = Object(reader.parseFile("sphere.obj"), Material(ReflectionType::diffuse, ColorD(1.0, 0, 0), ColorD(), 1.0, 0.0));
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

  obj = Object(reader.parseFile("sphere.obj"), Material(ReflectionType::specular, ColorD(), ColorD(), 0.5, 0.5));
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

  // right
  obj = Object(reader.parseFile("cube.obj"), Material(ReflectionType::diffuse, ColorD(1.0, 0, 0), ColorD(), 1.0, 0.0));
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
  obj = Object(reader.parseFile("cube.obj"), Material(ReflectionType::diffuse, ColorD(0, 0, 1.0), ColorD(), 1.0, 0.0));
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
  obj = Object(reader.parseFile("cube.obj"), Material(ReflectionType::diffuse, ColorD(1.0, 1.0, 1.0), ColorD(), 1.0, 0.0));
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
  obj = Object(reader.parseFile("cube.obj"), Material(ReflectionType::diffuse, ColorD(0.0, 1.0, 0.0), ColorD(), 1.0, 0.0));
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
  obj = Object(reader.parseFile("cube.obj"), Material(ReflectionType::diffuse, ColorD(1.0, 1.0, 0.0), ColorD(), 1.0, 0.0));
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

  //camera = Camera(Vector3D(0, 15.0, -10.0), Vector3D::Normalise(Vector3D(0, -0.5, 1)), Vector3D(0, 1, 0));
  camera = Camera(Vector3D(0, 0, 4), Vector3D::Normalise(Vector3D(0, 0, -1)), Vector3D(0, 1, 0));
}