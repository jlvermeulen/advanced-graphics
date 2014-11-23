//---------------------------------------------------------------------------------------
inline Vector3D Vector3D::operator+(const Vector3D& rhs)
{
	Vector3D result(*this);

	result.X += rhs.X;
	result.Y += rhs.Y;
	result.Z += rhs.Z;

	return result;
}

//--------------------------------------------------------------------------------------
inline Vector3D Vector3D::operator+(const double& rhs)
{
	Vector3D result(*this);

	result.X += rhs;
	result.Y += rhs;
	result.Z += rhs;

	return result;
}

//--------------------------------------------------------------------------------------
inline Vector3D operator+(const double& lhs, const Vector3D& rhs)
{
	Vector3D result(rhs);

	result.X += lhs;
	result.Y += lhs;
	result.Z += lhs;

	return result;
}

//--------------------------------------------------------------------------------------
inline Vector3D Vector3D::operator-()
{
	Vector3D result(*this);

	result.X *= -1;
	result.Y *= -1;
	result.Z *= -1;

	return result;
}

//--------------------------------------------------------------------------------------
inline Vector3D Vector3D::operator-(const Vector3D& rhs)
{
	Vector3D result(*this);

	result.X -= rhs.X;
	result.Y -= rhs.Y;
	result.Z -= rhs.Z;

	return result;
}

//--------------------------------------------------------------------------------------
inline Vector3D Vector3D::operator-(const double& rhs)
{
	Vector3D result(*this);

	result.X -= rhs;
	result.Y -= rhs;
	result.Z -= rhs;

	return result;
}

//--------------------------------------------------------------------------------------
inline Vector3D operator-(const double& lhs, const Vector3D& rhs)
{
	Vector3D result(rhs);

	result.X -= lhs;
	result.Y -= lhs;
	result.Z -= lhs;

	return result;
}

//--------------------------------------------------------------------------------------
inline Vector3D Vector3D::operator*(const double& rhs)
{
	Vector3D result(*this);

	result.X *= rhs;
	result.Y *= rhs;
	result.Z *= rhs;

	return result;
}

//--------------------------------------------------------------------------------------
inline Vector3D operator*(const double& lhs, const Vector3D& rhs)
{
	Vector3D result(rhs);

	result.X *= lhs;
	result.Y *= lhs;
	result.Z *= lhs;

	return result;
}

//--------------------------------------------------------------------------------------
inline Vector3D operator*(const Matrix3x3D& lhs, const Vector3D& rhs)
{
	Vector3D result(rhs);

	result.X = result.Dot(Vector3D(lhs[0]));
	result.Y = result.Dot(Vector3D(lhs[1]));
	result.Z = result.Dot(Vector3D(lhs[2]));

	return result;
}

//--------------------------------------------------------------------------------------
inline Vector3D Vector3D::operator/(const double& rhs)
{
	Vector3D result(*this);

	result.X /= rhs;
	result.Y /= rhs;
	result.Z /= rhs;

	return result;
}

//--------------------------------------------------------------------------------------
inline void Vector3D::operator+=(const Vector3D& rhs)
{
	this->X += rhs.X;
	this->Y += rhs.Y;
	this->Z += rhs.Z;
}

//--------------------------------------------------------------------------------------
inline void Vector3D::operator+=(const double& rhs)
{
	this->X += rhs;
	this->Y += rhs;
	this->Z += rhs;
}

//--------------------------------------------------------------------------------------
inline void Vector3D::operator-=(const Vector3D& rhs)
{
	this->X -= rhs.X;
	this->Y -= rhs.Y;
	this->Z -= rhs.Z;
}

//--------------------------------------------------------------------------------------
inline void Vector3D::operator-=(const double& rhs)
{
	this->X -= rhs;
	this->Y -= rhs;
	this->Z -= rhs;
}

//--------------------------------------------------------------------------------------
inline void Vector3D::operator*=(const double& rhs)
{
	this->X *= rhs;
	this->Y *= rhs;
	this->Z *= rhs;
}

//--------------------------------------------------------------------------------------
inline void Vector3D::operator/=(const double& rhs)
{
	this->X /= rhs;
	this->Y /= rhs;
	this->Z /= rhs;
}

//--------------------------------------------------------------------------------------
inline bool Vector3D::operator==(const Vector3D& rhs)
{
	return this->X == rhs.X && this->Y == rhs.Y && this->Z == rhs.Z;
}

//--------------------------------------------------------------------------------------
inline bool Vector3D::operator!=(const Vector3D& rhs)
{
	return this->X != rhs.X || this->Y != rhs.Y || this->Z != rhs.Z;
}