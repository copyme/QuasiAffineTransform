#ifndef _PAVING_H_
#define	_PAVING_H_

#include <vector>

#include "vector2d.hpp"

class Paving
{
 private:
  std::vector<Vector2D> m_points;
	
 public:
	
  unsigned int size() const;
  void addPoint(const Vector2D);
	
  const Vector2D operator[](unsigned int) const;
  Vector2D& operator[](unsigned int);
	
};

#endif //_PAVING_H_
