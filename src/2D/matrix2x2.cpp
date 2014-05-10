/**
 * @file   matrix2x2.cpp
 * @author Valentin Blot <valentin@Valentin-E21ter>
 * @date   Thu Jun 19 13:12:30 2008
 * 
 * @brief  2-dimensional matrix
 * 
 * 
 */

#include "matrix2x2.hpp"

using namespace std;

/** 
 * default constructor
 * 
 */
Matrix2x2::Matrix2x2()
{
}

/** 
 * cnstruct from 2 vectors
 * 
 * @param c1 first column
 * @param c2 second column
 */
Matrix2x2::Matrix2x2(const Vector2D c1, const Vector2D c2)
{
  m_tab[0] = c1;
  m_tab[1] = c2;
}

/** 
 * construct from 4 integers
 * 
 * @param v1 upper-left coefficient
 * @param v2 upper-right coefficient
 * @param v3 down-left coefficient
 * @param v4 down-right coefficient
 */
Matrix2x2::Matrix2x2(int v1, int v2, int v3, int v4)
{
  m_tab[0][0] = v1;
  m_tab[1][0] = v2;
  m_tab[0][1] = v3;
  m_tab[1][1] = v4;
}

/** 
 * computes determinant of the matrix
 * 
 * 
 * @return the determinant
 */
int Matrix2x2::det() const
{
  return m_tab[0][0] * m_tab[1][1] - m_tab[1][0] * m_tab[0][1] ;
}

/** 
 * computes the opposite of a matrix
 * 
 * 
 * @return the opposite of the matrix
 */
const Matrix2x2 Matrix2x2::operator-() const
{
  return Matrix2x2(-m_tab[0], -m_tab[1]);
}

/** 
 * multiply a matrix and a value
 * 
 * @param value 
 * 
 * @return 
 */
const Matrix2x2 Matrix2x2::operator*(int value) const
{
  return Matrix2x2(m_tab[0][0] * value, m_tab[1][0] * value, m_tab[0][1] * value, m_tab[1][1] * value); 
}

/** 
 * multiply current matrix by a value
 * 
 * @param value the other matrix
 * 
 * @return the modified matrix
 */
Matrix2x2& Matrix2x2::operator*=(int value)
{
  *this = *this * value;
  return *this;
}

/** 
 * multiply a matrix and a vector
 * 
 * @param vect the vector
 * 
 * @return the vector result of the product
 */
const Vector2D Matrix2x2::operator*(const Vector2D vect) const
{
  return m_tab[0] * vect[0] + m_tab[1] * vect[1];
}

/** 
 * multiply two matrix
 * 
 * @param mat the other matrix
 * 
 * @return the product
 */
const Matrix2x2 Matrix2x2::operator*(const Matrix2x2 mat) const
{
  return Matrix2x2(*this * mat[0], *this * mat[1]);
}

/** 
 * multiply current matrix by another matrix
 * 
 * @param mat the other matrix
 * 
 * @return the modified matrix
 */
Matrix2x2& Matrix2x2::operator*=(const Matrix2x2 mat)
{
  *this = *this * mat;
  return *this;
}

/** 
 * get column - read-only version
 * 
 * @param col column to get
 * 
 * @return column
 */
const Vector2D Matrix2x2::operator[](unsigned int col) const
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
Vector2D& Matrix2x2::operator[](unsigned int col)
{
  return m_tab[col];
}

void Matrix2x2::print()
{
  cout <<"{ "<< m_tab[0][0] <<" "<<m_tab[0][1] << "\n "<<"  "<<m_tab[1][0]<<" "<<m_tab[1][1]<<"}"<<endl;
}
