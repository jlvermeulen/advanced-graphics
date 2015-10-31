#pragma once

#include "Camera.h"
#include "ColorD.h"
#include "Material.h"
#include "Object.h"
#include "Octree.h"
#include "Triangle.h"

#include <vector>
#include <utility>
#include <tuple>
#include <random>

typedef unsigned char uchar;

struct Intersection
{
	Intersection(const Ray& ray, double time, const Triangle& hit, const Material& material) :
		hitPoint(ray.Origin + time * ray.Direction),
		hit(hit),
		hitMaterial(material)
	{
	}

	Triangle hit;
	Material hitMaterial;
	Vector3D hitPoint;
};

struct Light
{
	Light(const Vector3D& position, const ColorD& color) :
		position(position),
		color(color)
	{
	}

	Vector3D position;
	ColorD color;
};

class Scene
{
public:
	Scene();
	~Scene();

	void PreRender(int minTriangles, int maxDepth);
	void Render(unsigned char* imageData, int samplesPerPixel, double sigma, bool useDoF);
	void PostRender();

	void LoadDefaultScene();
	void LoadDefaultScene2();
	void LoadDefaultScene3();
	void LoadDefaultScene4();
	void LoadDefaultScene5();
	void LoadDefaultScene6();
	void Clear();

private:
	ColorD TraceRay(const Ray& ray, bool nee);
	ColorD ComputeRadiance(const Vector3D& point, const Vector3D& in, Triangle* triangle, const Material& material, unsigned int depth, bool nee);
	ColorD DirectIllumination(const Vector3D& point, const Vector3D& in, const Vector3D& normal, const Material& material);
	ColorD IndirectIllumination(Vector3D point, const Vector3D& in, const Vector3D& normal, const Material& material, unsigned int depth, bool nee);

	void TracePixels(std::pair<ColorD, double>* pixelData, int samplesPerPixel, double sigma, bool useDoF);

	double GaussianWeight(double dx, double dy, double sigma) const;

	std::pair<Ray, Triangle*>* SampleLight(const Vector3D& hitPoint, double& weight);
	Triangle* FirstHitInfo(const Ray& ray, double& time, Material& mat) const;

public:
	Camera camera;
	std::vector<Object*> objects;
	std::vector<Object*> lights;
	std::uniform_real_distribution<double> dist;
	std::mt19937 gen;
};