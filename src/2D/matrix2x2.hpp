#ifndef _MATRIX2X2_H_
#define _MATRIX2X2_H_

#include "vector2d.hpp"
#include <iostream>

class Matrix2x2
{
 private:
  Vector2D m_tab[2];
 public:
  Matrix2x2();
  Matrix2x2(int, int, int, int);
  Matrix2x2(const Vector2D, const Vector2D);
	
  int det() const;
	
  const Matrix2x2 operator-() const;
  const Matrix2x2 operator*(int) const;
  Matrix2x2& operator*=(int);
  const Matrix2x2 operator*(const Matrix2x2) const;
  const Vector2D operator*(const Vector2D) const;
  Matrix2x2& operator*=(const Matrix2x2);

  const Vector2D operator[](unsigned int) const;
  Vector2D& operator[](unsigned int);

   
   void print();
   
};


#endif
