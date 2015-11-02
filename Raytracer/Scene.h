#pragma once

#include "Camera.h"
#include "Color3F.h"
#include "Material.h"
#include "Object.h"
#include "Octree.h"
#include "Triangle.h"

#include <vector>
#include <utility>
#include <tuple>
#include <random>
#include <stdexcept>

typedef unsigned char uchar;

struct Intersection
{
	Intersection(const Ray& ray, float time, const Triangle& hit, const Material& material) :
		hitPoint(ray.Origin + time * ray.Direction),
		hit(hit),
		hitMaterial(material)
	{
	}

	Triangle hit;
	Material hitMaterial;
	Vector3F hitPoint;
};

struct Light
{
	Light(const Vector3F& position, const Color3F& color) :
		position(position),
		color(color)
	{
	}

	Vector3F position;
	Color3F color;
};

class Scene
{
public:
	Scene();
	~Scene();

	void PreRender();
	void Render(unsigned char* imageData, int samplesPerPixel, float sigma, bool useDoF);
	void PostRender();

	void LoadDefaultScene();
	void LoadDefaultScene2();
	void LoadDefaultScene3();
	void LoadDefaultScene4();
	void LoadDefaultScene5();
	void LoadDefaultScene6();
	void Clear();

private:
	Color3F TraceRay(const Ray& ray, bool nee);
	Color3F ComputeRadiance(const Vector3F& point, const Vector3F& in, Triangle* triangle, const Material& material, unsigned int depth, bool nee);
	Color3F DirectIllumination(const Vector3F& point, const Vector3F& in, const Vector3F& normal, const Material& material);
	Color3F IndirectIllumination(Vector3F point, const Vector3F& in, const Vector3F& normal, const Material& material, unsigned int depth, bool nee);

	void TracePixels(std::pair<Color3F, float>* pixelData, int samplesPerPixel, float sigma, bool useDoF);

	float GaussianWeight(float dx, float dy, float sigma) const;

	std::pair<Ray, Triangle*>* SampleLight(const Vector3F& hitPoint, float& weight);
	Triangle* FirstHitInfo(const Ray& ray, float& time, Material& mat) const;

public:
	Camera camera;
	std::vector<Object*> objects;
	std::vector<Object*> lights;
	unsigned int lightCount;
	std::uniform_real_distribution<float> dist;
	std::mt19937 gen;
};