/**
 * @file   mathematic.cpp
 * @author Valentin Blot <valentin@Valentin-E21ter>
 * @date   Thu Jun 19 14:01:00 2008
 * 
 * @brief  mathematic functions
 * 
 * 
 */

#include "mathematic.hpp"

namespace Math
{

  /** 
   * computes the greatest common divider of 2 integers
   * 
   * @param a first integer
   * @param b second integer
   * 
   * @return 
   */
  int pgcd( int a, int b )
  {
    a = abs( a );
    b = abs( b );
    int t;
    while( b != 0 )
      {
	t = b ;
	b = a % b ;
	a = t ;
      }
    return a;
  }

  /** 
   * computes v0 and v1 such as v0 * a + v1 * b = gcd(a, b)
   * 
   * @param a the first integer
   * @param b the second integer
   * 
   * @return the vector (v0, v1)
   */
  couple extendedEuclidean(int a, int b)
  {
    couple v;

    if(b == 0)
      {
	v[0] = a < 0? -1: 1;
	v[1] = 0;
      }
    else
      {
	couple vt = extendedEuclidean(b, a % b);
	v[0] = vt[1];
	v[1] = vt[0] - a / b * vt[1];
      }
    return v;
  }

  /** 
   * quotient of the Euclidean division of 2 integers with positive remain
   * 
   * @param a first integer
   * @param b second integer
   * 
   * @return 
   */
  int intDiv( int a, int b )
  {
    int q = a / b;
    if(a - b * q < 0)
      return q - 1;
    return q;
  }

  /** 
   * remain of the Euclidean division of 2 integers with positive remain
   * 
   * @param a first integer
   * @param b second integer
   * 
   * @return 
   */
  int mod( int a, int b )
  {
    return a - b * intDiv( a, b ) ;
  }

} //namespace Math
