#ifndef _VECTOR2D_H_
#define _VECTOR2D_H_

#include "../common/mathematic.hpp"

class Vector2D
{
 private:
  int m_tab[2];

 public:
  Vector2D();
  Vector2D(int, int);

  int operator[](unsigned int) const;
  int &operator[](unsigned int);
  const Vector2D operator+(const Vector2D) const;
  const Vector2D operator-(const Vector2D) const;
  const Vector2D operator*(int) const;
  const Vector2D operator/(int) const;
  const Vector2D operator%(int) const;
  const Vector2D operator-() const;
  Vector2D& operator+=(const Vector2D);
  Vector2D& operator-=(const Vector2D);

};

#endif // _VECTOR2D_H_
