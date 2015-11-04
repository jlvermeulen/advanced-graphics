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
#define MAXLIGHTS 32
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

	union{ float  posv0X[MAXLIGHTS]; __m256 posv0X8[MAXLIGHTS / NROFLANES]; };
	union{ float  posv1X[MAXLIGHTS]; __m256 posv1X8[MAXLIGHTS / NROFLANES]; };
	union{ float  posv2X[MAXLIGHTS]; __m256 posv2X8[MAXLIGHTS / NROFLANES]; };
	union{ float  posv0Y[MAXLIGHTS]; __m256 posv0Y8[MAXLIGHTS / NROFLANES]; };
	union{ float  posv1Y[MAXLIGHTS]; __m256 posv1Y8[MAXLIGHTS / NROFLANES]; };
	union{ float  posv2Y[MAXLIGHTS]; __m256 posv2Y8[MAXLIGHTS / NROFLANES]; };
	union{ float  posv0Z[MAXLIGHTS]; __m256 posv0Z8[MAXLIGHTS / NROFLANES]; };
	union{ float  posv1Z[MAXLIGHTS]; __m256 posv1Z8[MAXLIGHTS / NROFLANES]; };
	union{ float  posv2Z[MAXLIGHTS]; __m256 posv2Z8[MAXLIGHTS / NROFLANES]; };

	union{ float  norv0X[MAXLIGHTS]; __m256 norv0X8[MAXLIGHTS / NROFLANES]; };
	union{ float  norv1X[MAXLIGHTS]; __m256 norv1X8[MAXLIGHTS / NROFLANES]; };
	union{ float  norv2X[MAXLIGHTS]; __m256 norv2X8[MAXLIGHTS / NROFLANES]; };
	union{ float  norv0Y[MAXLIGHTS]; __m256 norv0Y8[MAXLIGHTS / NROFLANES]; };
	union{ float  norv1Y[MAXLIGHTS]; __m256 norv1Y8[MAXLIGHTS / NROFLANES]; };
	union{ float  norv2Y[MAXLIGHTS]; __m256 norv2Y8[MAXLIGHTS / NROFLANES]; };
	union{ float  norv0Z[MAXLIGHTS]; __m256 norv0Z8[MAXLIGHTS / NROFLANES]; };
	union{ float  norv1Z[MAXLIGHTS]; __m256 norv1Z8[MAXLIGHTS / NROFLANES]; };
	union{ float  norv2Z[MAXLIGHTS]; __m256 norv2Z8[MAXLIGHTS / NROFLANES]; };

	union{ float  prev0X[MAXLIGHTS]; __m256 prev0X8[MAXLIGHTS / NROFLANES]; };
	union{ float  prev0Y[MAXLIGHTS]; __m256 prev0Y8[MAXLIGHTS / NROFLANES]; };
	union{ float  prev0Z[MAXLIGHTS]; __m256 prev0Z8[MAXLIGHTS / NROFLANES]; };
	union{ float  prev1X[MAXLIGHTS]; __m256 prev1X8[MAXLIGHTS / NROFLANES]; };
	union{ float  prev1Y[MAXLIGHTS]; __m256 prev1Y8[MAXLIGHTS / NROFLANES]; };
	union{ float  prev1Z[MAXLIGHTS]; __m256 prev1Z8[MAXLIGHTS / NROFLANES]; };

	union{ float  area[MAXLIGHTS]; __m256 area8[MAXLIGHTS / NROFLANES]; };
	union{ float  emissionR[MAXLIGHTS]; __m256 emissionR8[MAXLIGHTS / NROFLANES]; };
	union{ float  invDenom[MAXLIGHTS]; __m256 invDenom8[MAXLIGHTS / NROFLANES]; };
	union{ float  d00[MAXLIGHTS]; __m256 d008[MAXLIGHTS / NROFLANES]; };
	union{ float  d01[MAXLIGHTS]; __m256 d018[MAXLIGHTS / NROFLANES]; };
	union{ float  d11[MAXLIGHTS]; __m256 d118[MAXLIGHTS / NROFLANES]; };
	Triangle* lightTriangles[MAXLIGHTS];

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