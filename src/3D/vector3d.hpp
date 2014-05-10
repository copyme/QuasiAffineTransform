#ifndef _VECTOR3D_H_
#define _VECTOR3D_H_

#include "../common/mathematic.hpp"

class Vector3D
{
 private:
  int m_tab[3];

 public:
  Vector3D();
  Vector3D(int, int, int);

  int operator[](unsigned int) const;
  int &operator[](unsigned int);
  const Vector3D operator+(const Vector3D) const;
  const Vector3D operator-(const Vector3D) const;
  const Vector3D operator*(int) const;
  const Vector3D operator/(int) const;
  Vector3D& operator+=(const Vector3D);
  Vector3D& operator-=(const Vector3D);
  const Vector3D operator-() const;

};

#endif // _VECTOR3D_H_
