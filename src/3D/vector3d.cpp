/**
 * @file   vector2d.cpp
 * @author Valentin Blot <valentin@Valentin-E21ter>
 * @date   Thu Jun 19 11:22:30 2008
 * 
 * @brief  2-dimensinal vectors
 * 
 * 
 */
#include "vector3d.hpp"

/** 
 * default constructor
 * 
 */
Vector3D::Vector3D()
{
}

/** 
 * construct from 3 integers
 * 
 * @param v1 first coordinate
 * @param v2 second coordinate
 * @param v3 third coordinate
 */
Vector3D::Vector3D(int v1, int v2, int v3)
{
  m_tab[0] = v1;
  m_tab[1] = v2;
  m_tab[2] = v3;
}

/** 
 * get coordinate - read-only version
 * 
 * @param position coordinate to get
 * 
 * @return coordinate
 */
int Vector3D::operator[](unsigned int position) const
{
  return m_tab[position];
}

/** 
 * get coordinate - read-write version
 * 
 * @param position coordinate to get
 * 
 * @return coordinate
 */
int &Vector3D::operator[](unsigned int position)
{
  return m_tab[position] ;
}

/** 
 * add two vectors
 * 
 * @param vect the other vector
 * 
 * @return the sum
 */
const Vector3D Vector3D::operator+(const Vector3D vect) const
{
  return Vector3D(m_tab[0] + vect.m_tab[0], m_tab[1] + vect.m_tab[1], m_tab[2] + vect.m_tab[2]);
}

/** 
 * substract two vectors
 * 
 * @param vect the other vector
 * 
 * @return the substraction
 */
const Vector3D Vector3D::operator-(const Vector3D vect) const
{
  return *this + (-vect);
}


/** 
 * multiply a vector and a value
 * 
 * @param value the value
 * 
 * @return the product
 */
const Vector3D Vector3D::operator*(int value) const
{
  return Vector3D(m_tab[0] * value, m_tab[1] * value, m_tab[2] * value);
}

/** 
 * divides a vector by a value
 * 
 * @param value the value
 * 
 * @return the result
 */
const Vector3D Vector3D::operator/(int value) const
{
  return Vector3D(Math::intDiv(m_tab[0], value),
		  Math::intDiv(m_tab[1], value),
		  Math::intDiv(m_tab[2], value));
}

/** 
 * add another vector to current vector
 * 
 * @param vect the other vector
 * 
 * @return the modified vector
 */
Vector3D &Vector3D::operator+=(const Vector3D vect)
{
  *this = *this + vect;
  return *this;
}


/** 
 * substract another vector to current vector
 * 
 * @param vect the other vector
 * 
 * @return the modified vector
 */
Vector3D& Vector3D::operator-=(const Vector3D vect)
{
  *this = *this - vect;
  return *this;
}

/** 
 * computes the opposite of a vector
 * 
 * 
 * @return the opposite of the vector
 */
const Vector3D Vector3D::operator-() const
{
  return Vector3D(-m_tab[0], -m_tab[1], -m_tab[2]);
}
