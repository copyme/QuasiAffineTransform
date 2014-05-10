/**
 * @file   vector2d.cpp
 * @author Valentin Blot <valentin@Valentin-E21ter>
 * @date   Thu Jun 19 11:22:30 2008
 * 
 * @brief  2-dimensinal vectors
 * 
 * 
 */
#include "vector2d.hpp"

/** 
 * default constructor
 * 
 */
Vector2D::Vector2D()
{
}

/** 
 * construct from 2 integers
 * 
 * @param v1 first coordinate
 * @param v2 second coordinate
 */
Vector2D::Vector2D(int v1, int v2)
{
  m_tab[0] = v1;
  m_tab[1] = v2;
}

/** 
 * get coordinate - read-only version
 * 
 * @param position coordinate to get
 * 
 * @return coordinate
 */
int Vector2D::operator[](unsigned int position) const
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
int &Vector2D::operator[](unsigned int position)
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
const Vector2D Vector2D::operator+(const Vector2D vect) const
{
  return Vector2D(m_tab[0] + vect.m_tab[0], m_tab[1] + vect.m_tab[1]);
}

/** 
 * substract two vectors
 * 
 * @param vect the other vector
 * 
 * @return the substraction
 */
const Vector2D Vector2D::operator-(const Vector2D vect) const
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
const Vector2D Vector2D::operator*(int value) const
{
  return Vector2D(m_tab[0] * value, m_tab[1] * value);
}

/** 
 * divides a vector by a value
 * 
 * @param value the value
 * 
 * @return the result
 */
const Vector2D Vector2D::operator/(int value) const
{
  return Vector2D(Math::intDiv(m_tab[0], value), Math::intDiv(m_tab[1], value));
}

/** 
* Return the remainder of the division of a vector by a value
* 
* @param value the value
* 
* @return the result
*/
const Vector2D Vector2D::operator%(int value) const
{
  return Vector2D(Math::mod(m_tab[0], value), Math::mod(m_tab[1], value));
}

/** 
 * computes the opposite of a vector
 * 
 * 
 * @return the opposite vector
 */
const Vector2D Vector2D::operator-() const
{
  return Vector2D(-m_tab[0], -m_tab[1]);
}

/** 
 * add another vector to current vector
 * 
 * @param vect the other vector
 * 
 * @return the modified vector
 */
Vector2D &Vector2D::operator+=(const Vector2D vect)
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
Vector2D &Vector2D::operator-=(const Vector2D vect)
{
  *this = *this - vect;
  return *this;
}
