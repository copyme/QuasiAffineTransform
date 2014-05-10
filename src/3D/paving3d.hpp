#ifndef _PAVING3D_H_
#define	_PAVING3D_H_

#include <vector>

#include "vector3d.hpp"

class Paving3D
{
 private:
  std::vector<Vector3D> m_points;
	
 public:
	
  unsigned int size() const;
  void addPoint(const Vector3D);

  const Vector3D operator[](unsigned int) const;
  Vector3D& operator[](unsigned int);
	
};

#endif //_PAVING3D_H_
