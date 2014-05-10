#ifndef _MATHEMATIC_H_
#define	_MATHEMATIC_H_

#include <cstdarg>
#include <cstdlib>

namespace Math
{
  struct couple
  {
    int a;
    int b;
    int &operator[](int i)
    {
      if(i == 0)
	return a;
      return b;
    }
  };

  int pgcd(int, int);

  couple extendedEuclidean(int, int);

  int intDiv(int, int);
  int mod(int, int);


} //namespace Math

#endif //_MATHEMATIC_H_
