#include "BBox.h"
#include <algorithm>
#include "BVHRay.h"
#include <iostream>
#include <fstream>

BBox::BBox(const Vector3& min, const Vector3& max)
  : min(min), max(max) { extent = max - min; }

BBox::BBox(const Vector3& p)
  : min(p), max(p) { extent = max - min; }

  void BBox::expandToInclude(const Vector3& p) {
    min = ::min(min, p);
    max = ::max(max, p);
    extent = max - min;
  }

void BBox::expandToInclude(const BBox& b) {
  min = ::min(min, b.min);
  max = ::max(max, b.max);
  extent = max - min;
}

uint32_t BBox::maxDimension() const {
  uint32_t result = 0;
  if(extent.y > extent.x) result = 1;
  if(extent.z > extent.y) result = 2;
  return result;
}

float BBox::surfaceArea() const {
  return 2.f*( extent.x*extent.z + extent.x*extent.y + extent.y*extent.z );
}

// http://www.flipcode.com/archives/SSE_RayBox_Intersection_Test.shtml
// turn those verbose intrinsics into something readable.
#define loadps(mem)		_mm_loadu_ps((const float * const)(mem))
#define storess(ss,mem)		_mm_store_ss((float * const)(mem),(ss))
#define minss			_mm_min_ss
#define maxss			_mm_max_ss
#define minps			_mm_min_ps
#define maxps			_mm_max_ps
#define mulps			_mm_mul_ps
#define subps			_mm_sub_ps
#define rotatelps(ps)		_mm_shuffle_ps((ps),(ps), 0x39)	// a,b,c,d -> b,c,d,a
#define muxhps(low,high)	_mm_movehl_ps((low),(high))	// low{a,b,c,d}|high{e,f,g,h} = {c,d,g,h}
static const float flt_plus_inf = -logf(0);	// let's keep C and C++ compilers happy.
static const float 
  ps_cst_plus_inf[4] = {  flt_plus_inf,  flt_plus_inf,  flt_plus_inf,  flt_plus_inf },
  ps_cst_minus_inf[4] = { -flt_plus_inf, -flt_plus_inf, -flt_plus_inf, -flt_plus_inf };
bool BBox::intersect(const Ray& ray, float *tnear, float *tfar) const {

	std::cout << "start";
		Vector3D o = ray.Origin;
	Vector3D d = ray.Direction;
	__declspec(align(32))Vector3 oo(o.X, o.Y, o.Z);
	__declspec(align(32))Vector3 dd(d.X, d.Y, d.Z);
	__declspec(align(32))BVHRay bvhray(oo, dd);
	
  // you may already have those values hanging around somewhere
	__declspec(align(32)) const __m128
    plus_inf	= loadps(ps_cst_plus_inf),
              minus_inf	= loadps(ps_cst_minus_inf);

	std::cout << bvhray.o.x;
  // use whatever's apropriate to load.
  __declspec(align(32)) const __m128 box_min = min.m128;//loadps(&min);
  __declspec(align(32)) const __m128 box_max = max.m128;// loadps(&max);
  __declspec(align(32))const __m128 pos = bvhray.o.m128;// loadps(&bvhray.o);
  __declspec(align(32)) const __m128 inv_dir = bvhray.inv_d.m128;// loadps(&bvhray.inv_d);

  std::cout << pos.m128_f32;
  std::cout << box_min.m128_f32;
  std::cout << inv_dir.m128_f32;
  // use a div if inverted directions aren't available
  __declspec(align(32))const __m128 h1 = subps(box_min, pos);
  __declspec(align(32))const __m128 l1 = mulps(h1, inv_dir);
  __declspec(align(32))const __m128 l2 = mulps(subps(box_max, pos), inv_dir);

  // the order we use for those min/max is vital to filter out
  // NaNs that happens when an inv_dir is +/- inf and
  // (box_min - pos) is 0. inf * 0 = NaN
  const __m128 filtered_l1a = minps(l1, plus_inf);
  const __m128 filtered_l2a = minps(l2, plus_inf);

  const __m128 filtered_l1b = maxps(l1, minus_inf);
  const __m128 filtered_l2b = maxps(l2, minus_inf);

  // now that we're back on our feet, test those slabs.
  __m128 lmax = maxps(filtered_l1a, filtered_l2a);
  __m128 lmin = minps(filtered_l1b, filtered_l2b);

  // unfold back. try to hide the latency of the shufps & co.
  const __m128 lmax0 = rotatelps(lmax);
  const __m128 lmin0 = rotatelps(lmin);
  lmax = minss(lmax, lmax0);
  lmin = maxss(lmin, lmin0);

  const __m128 lmax1 = muxhps(lmax,lmax);
  const __m128 lmin1 = muxhps(lmin,lmin);
  lmax = minss(lmax, lmax1);
  lmin = maxss(lmin, lmin1);

  const bool ret = _mm_comige_ss(lmax, _mm_setzero_ps()) & _mm_comige_ss(lmax,lmin);

  storess(lmin, tnear);
  storess(lmax, tfar);

  return  ret;
}
