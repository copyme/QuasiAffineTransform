#ifndef _MATRIX3X3_H_
#define _MATRIX3X3_H_

#include "vector3d.hpp"

class Matrix3x3
{
 private:
  Vector3D m_tab[3];
 public:
  Matrix3x3();
  Matrix3x3(int, int, int, int, int, int, int, int, int);
  Matrix3x3(const Vector3D, const Vector3D, const Vector3D);
	
  int det() const;
	
  const Matrix3x3 operator*(int) const;
  Matrix3x3 &operator*=(int);
  const Matrix3x3 operator*(const Matrix3x3) const;
  Matrix3x3 &operator*=(const Matrix3x3);
  const Vector3D operator*(const Vector3D) const;
  const Matrix3x3 operator-() const;

  const Vector3D operator[](unsigned int) const;
  Vector3D& operator[](unsigned int);

};


/** 
 * default constructor
 * 
 */
inline Matrix3x3::Matrix3x3()
{
}

/** 
 * construct from 9 integers
 * 
 * @param v1 upper-left coefficient
 * @param v2 upper coefficient
 * @param v3 upper-right coefficient
 * @param v4 left coefficient
 * @param v5 center coefficient
 * @param v6 right coefficient
 * @param v7 down-left coefficient
 * @param v8 down coefficient
 * @param v9 down-right coefficient
 */
inline Matrix3x3::Matrix3x3(int v1, int v2, int v3, int v4, int v5, int v6, int v7, int v8, int v9)
{
  m_tab[0][0] = v1;
  m_tab[1][0] = v2;
  m_tab[2][0] = v3;
  m_tab[0][1] = v4;
  m_tab[1][1] = v5;
  m_tab[2][1] = v6;
  m_tab[0][2] = v7;
  m_tab[1][2] = v8;
  m_tab[2][2] = v9;
}

/** 
 * construct from 3 vectors
 * 
 * @param v1 first column
 * @param v2 second column
 * @param v3 third column
 */
inline Matrix3x3::Matrix3x3(const Vector3D v1, const Vector3D v2, const Vector3D v3)
{
  m_tab[0] = v1;
  m_tab[1] = v2;
  m_tab[2] = v3;
}

/** 
 * computes determinant of the matrix
 * 
 * 
 * @return the determinant
 */
inline int Matrix3x3::det() const
{
  return m_tab[0][0] * m_tab[1][1] * m_tab[2][2] + m_tab[1][0] * m_tab[2][1] * m_tab[0][2] + m_tab[2][0] * m_tab[0][1] * m_tab[1][2]
    - m_tab[2][0] * m_tab[1][1] * m_tab[0][2] - m_tab[1][0] * m_tab[0][1] * m_tab[2][2] - m_tab[0][0] * m_tab[2][1] * m_tab[1][2];
}

/** 
 * multiply a matrix and a value
 * 
 * @param value the value
 * 
 * @return 
 */
inline const Matrix3x3 Matrix3x3::operator*(int value) const
{
  return Matrix3x3(m_tab[0] * value, m_tab[1] * value, m_tab[2] * value);
}

/** 
 * multiply current matrix by a value
 * 
 * @param value the value
 * 
 * @return the modified matrix
 */
inline Matrix3x3 &Matrix3x3::operator*=(int value)
{
  *this = *this * value;
  return *this;
}

/** 
 * multiply two matrix
 * 
 * @param value the other matrix
 * 
 * @return the product
 */
inline const Matrix3x3 Matrix3x3::operator*(const Matrix3x3 mat) const
{
  return Matrix3x3(m_tab[0] * mat.m_tab[0][0] + m_tab[1] * mat.m_tab[0][1] + m_tab[2] * mat.m_tab[0][2],
		   m_tab[0] * mat.m_tab[1][0] + m_tab[1] * mat.m_tab[1][1] + m_tab[2] * mat.m_tab[1][2],
		   m_tab[0] * mat.m_tab[2][0] + m_tab[1] * mat.m_tab[2][1] + m_tab[2] * mat.m_tab[2][2]);
}

/** 
 * multiply current matrix by another matrix
 * 
 * @param value the other matrix
 * 
 * @return the modified matrix
 */
inline Matrix3x3 &Matrix3x3::operator*=(const Matrix3x3 mat)
{
  *this = *this * mat;
  return *this;
}

/** 
 * multiply a matrix and a vector
 * 
 * @param vect the vector
 * 
 * @return the vector result of the product
 */
inline const Vector3D Matrix3x3::operator*(const Vector3D vect) const
{
  return m_tab[0] * vect[0] + m_tab[1] * vect[1] + m_tab[2] * vect[2];
}

/** 
 * computes the opposite of a matrix
 * 
 * 
 * @return the opposite of the matrix
 */
inline const Matrix3x3 Matrix3x3::operator-() const
{
  return Matrix3x3(-m_tab[0], -m_tab[1], -m_tab[2]);
}

/** 
 * get column - read-only version
 * 
 * @param col column to get
 * 
 * @return column
 */
inline const Vector3D Matrix3x3::operator[](unsigned int col) const
{
  return m_tab[col];
}

/** 
 * get column - read-write version
 * 
 * @param col column to get
 * 
 * @return column
 */
inline Vector3D& Matrix3x3::operator[](unsigned int col)
{
  return m_tab[col];
}


#endif
