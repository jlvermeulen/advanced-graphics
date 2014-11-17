#pragma once

struct Vector3D
{
public:
	double X, Y, Z;
	Vector3D() { }
	Vector3D(double x, double y, double z) : X(x), Y(y), Z(z) { }
	~Vector3D() { }

	inline double Dot(const Vector3D& other) { return Dot(*this, other); }
	inline Vector3D Cross(const Vector3D& other) { return Cross(*this, other); }
	inline double Length() { return sqrt(Dot(*this, *this)); }

	inline static double Dot(const Vector3D& lhs, const Vector3D& rhs) { return lhs.X * rhs.X + lhs.Y * rhs.Y + lhs.Z * rhs.Z; }
	inline static Vector3D Cross(const Vector3D& lhs, const Vector3D& rhs) { return Vector3D(lhs.Y * rhs.Z - lhs.Z * rhs.Y, lhs.Z * rhs.X - lhs.X * rhs.Z, lhs.X * rhs.Y - lhs.Y * rhs.X); }

	Vector3D& operator+=(const Vector3D& rhs)
	{
		X += rhs.X;
		Y += rhs.Y;
		Z += rhs.Z;
		return *this;
	}

	Vector3D& operator+=(const double& rhs)
	{
		X += rhs;
		Y += rhs;
		Z += rhs;
		return *this;
	}

	Vector3D& operator-=(const Vector3D& rhs)
	{
		X -= rhs.X;
		Y -= rhs.Y;
		Z -= rhs.Z;
		return *this;
	}

	Vector3D& operator-=(const double& rhs)
	{
		X -= rhs;
		Y -= rhs;
		Z -= rhs;
		return *this;
	}
	
	Vector3D& operator*=(const double& rhs)
	{
		X *= rhs;
		Y *= rhs;
		Z *= rhs;
		return *this;
	}

	Vector3D& operator/=(const double& rhs)
	{
		X /= rhs;
		Y /= rhs;
		Z /= rhs;
		return *this;
	}

private:
};

inline Vector3D& operator+(Vector3D lhs, const Vector3D& rhs) { lhs += rhs; return lhs; }
inline Vector3D& operator+(Vector3D lhs, const double& rhs) { lhs += rhs; return lhs; }
inline Vector3D& operator+(const double& lhs, Vector3D rhs) { rhs += lhs; return rhs; }

inline Vector3D& operator-(Vector3D lhs, const Vector3D& rhs) { lhs -= rhs; return lhs; }
inline Vector3D& operator-(Vector3D lhs, const double& rhs) { lhs -= rhs; return lhs; }
inline Vector3D& operator-(const double& lhs, Vector3D rhs) { rhs -= lhs; return rhs; }

inline Vector3D& operator*(Vector3D lhs, const double& rhs) { lhs *= rhs; return lhs; }
inline Vector3D& operator*(const double& lhs, Vector3D rhs) { rhs *= lhs; return rhs; }

inline Vector3D& operator/(Vector3D lhs, const double& rhs) { lhs /= rhs; return lhs; }

inline bool operator==(const Vector3D& lhs, const Vector3D& rhs) { return lhs.X == rhs.X && lhs.Y == rhs.Y && rhs.Z == lhs.Z; }
inline bool operator!=(const Vector3D& lhs, const Vector3D& rhs) { return !operator==(lhs, rhs); }