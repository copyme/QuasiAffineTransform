/**
 * @file   qat.cpp
 * @author Valentin Blot <valentin@Valentin-E21ter>
 * @date   Thu Jun 19 14:28:51 2008
 *
 * @brief  quasi-affine transforms utilities
 *
 *
 */

#include "qat.hpp"

/**
 * default constructor
 *
 */
QAT::QAT()
{
}

/**
 * construct from a matrix, an integer, and a vector
 *
 * @param matrix the matrix of the qat
 * @param omega the omega coefficient
 * @param vect the translation vector of the qat
 */
QAT::QAT ( Matrix2x2 matrix, int omega, Vector2D vect )
{
  m_matrix = matrix;
  m_omega = omega;
  m_vector = vect;
}

/**
 * determines if the qat is contracting
 *
 *
 * @return true if the qat is contracting
 */
bool QAT::isContracting() const
  {
    return m_omega * m_omega > m_matrix.det();
  }

/**
 * determines if the qat is inversible
 *
 *
 * @return true if the qat is inversible
 */
bool QAT::isInversible() const
  {
    int det = m_matrix.det();
    return det >= m_omega * std::max ( std::max ( abs ( m_matrix[1][1] ), abs ( m_matrix[0][1] ) ), abs ( m_matrix[1][0] ) + abs ( m_matrix[0][0] ) )
           || det >= m_omega * std::max ( std::max ( abs ( m_matrix[1][0] ), abs ( m_matrix[0][0] ) ), abs ( m_matrix[1][1] ) + abs ( m_matrix[0][1] ) );
  }

/**
 * determines the set of antecedents of a point by the qat - version using multiplications
 *
 * @param i X-axis coordinate of the point
 * @param j Y-axis coordinate of the point
 *
 * @return the paving containing the antecedents
 */
const Paving QAT::determinePaving ( int i, int j ) const
  {
    Paving paving;

    int A = -Math::intDiv ( -j * m_omega + f0, d1 );
    int B = -Math::intDiv ( - ( j + 1 ) * m_omega + f0, d1 );

    for ( int y = A; y < B; y++ )
      {
        int C = -Math::intDiv ( -i * m_omega + e0 + b1 * y, a1 );
        int D = -Math::intDiv ( - ( i + 1 ) * m_omega + e0 + b1 * y, a1 );
        for ( int x = C; x < D; x++ )
          paving.addPoint ( H * Vector2D ( x, y ) );
      }
    return paving;
  }

/**
 * determines the set of antecedents of a point by the qat - version using remains
 *
 * @param i X-axis coordinate of the point
 * @param j Y-axis coordinate of the point
 *
 * @return the paving containing the antecedents
 */
const Paving QAT::determinePavingNoMultiply ( int i, int j ) const
  {
    Paving paving;

    int A = -Math::intDiv ( -j * m_omega + f0, d1 );
    int B = -Math::intDiv ( - ( j + 1 ) * m_omega + f0, d1 );

    int C = -Math::intDiv ( -i * m_omega + e0 + b1 * A, a1 );
    int D = -Math::intDiv ( - ( i + 1 ) * m_omega + e0 + b1 * A, a1 );

    int Cp = Math::mod ( -i * m_omega + e0 + b1 * A, a1 );
    int Dp = Math::mod ( - ( i + 1 ) * m_omega + e0 + b1 * A, a1 );

    Vector2D X0 = H * Vector2D ( C, A );

    for ( int y = A; y < B; y++ )
      {
        Vector2D X = X0;
        for ( int x = C; x < D; x++ )
          {
            paving.addPoint ( X );
            X += H[0];
          }
        C -= b10;
        Cp += b11;
        X0 += H[1] - H0[0];
        if ( Cp >= a1 )
          {
            C--;
            Cp -= a1;
            X0 -= H[0];
          }
        D -= b10;
        Dp += b11;
        if ( Dp >= a1 )
          {
            D--;
            Dp -= a1;
          }
      }
    return paving;
  }

/**
 * determines the set of antecedents of a point by the qat - naive version
 *
 * @param i X-axis coordinate of the point
 * @param j Y-axis coordinate of the point
 *
 * @return the paving containing the antecedents
 */
const Paving QAT::determinePavingByBoundingRect ( int i, int j ) const
  {
    Paving paving;
    int delta = m_matrix.det();
    for ( int U = m_omega * i; U < m_omega * ( i + 1 ); U++ )
      for ( int V = m_omega * j; V < m_omega * ( j + 1 ); V++ )
        {
          int x = d0 * ( U - e0 ) - b0 * ( V - f0 );
          int y = -c0 * ( U - e0 ) + a0 * ( V - f0 );
          if ( x % delta == 0 && y % delta == 0 )
            {
              x /= delta;
              y /= delta;
              paving.addPoint ( Vector2D ( x, y ) );
            }
        }

    return paving;
  }

/**
 * adds a point to a paving in the set of pavings of the initial period
 *
 * @param i first coordinate of the paving
 * @param j second coordinate of the paving
 * @param P the point to add
 */
void QAT::setPaving ( const Vector2D I, const Vector2D P )
{
  if ( 0 <= I[0] && I[0] < alpha0 && 0 <= I[1] && I[1] < alpha1 )
    pavings[I[1] * alpha0 + I[0]].addPoint ( P );
}

/**
* adds remainders to a paving in the set of pavings of the initial period
*
* @param i first coordinate of the paving
* @param j second coordinate of the paving
* @param P the point to add
*/
void QAT::setPavingRemainder ( const Vector2D I, const Vector2D Rem, const Vector2D P )
{
  if ( 0 <= I[0] && I[0] < alpha0 && 0 <= I[1] && I[1] < alpha1 )
  {
      pavings[I[1] * alpha0 + I[0]].addPoint ( P );
      pavingsRemainder[I[1] * alpha0 + I[0]].addPoint ( Rem );
  }
}


/**
 * sets a paving in the set of pavings of the initial period
 *
 * @param i first coordinate of the paving
 * @param j second coordinate of the paving
 * @param P the paving
 */
void QAT::setPaving ( int i, int j, const Paving P )
{
  if ( 0 <= i && i < alpha0 && 0 <= j && j < alpha1 )
    pavings[j * alpha0 + i] = P;
}

/**
 * gets a paving in the set of pavings of the initial period
 *
 * @param i first coordinate of the paving
 * @param j second coordinate of the paving
 *
 * @return the paving
 */
const Paving QAT::getPaving ( int i, int j ) const
  {
    if ( 0 <= i && i < alpha0 && 0 <= j && j < alpha1 )
      return pavings[j * alpha0 + i];
    return Paving();
  }

/**
* gets a paving in the set of pavings of the initial period
*
* @param i first coordinate of the paving
* @param j second coordinate of the paving
*
* @return the paving
*/
const Paving QAT::getPavingRemainder ( int i, int j ) const
{
  if ( 0 <= i && i < alpha0 && 0 <= j && j < alpha1 )
    return pavingsRemainder[j * alpha0 + i];
  return Paving();
}



/**
* determines the set of pavings of the initial period - version with multiplications
*
*/
void QAT::determinePavings()
{
  int A = -Math::intDiv ( f0, d1 );
  for ( int y = A; y < A + m_omega * alpha1 / d1; y++ )
  {
    int C = -Math::intDiv ( e0 + b1 * y, a1 );
    for ( int x = C; x < C + m_omega * alpha0 / a1; x++ )
    {
      Vector2D X = H * Vector2D ( x, y );
      Vector2D I = calculate ( X );
      setPaving( I,  X );
    }
  }
}

/**
 * determines the set of pavings of the initial period - version with multiplications
 *
 */
void QAT::determinePavingsWithRemainders()
{
  int A = -Math::intDiv ( f0, d1 );
  for ( int y = A; y < A + m_omega * alpha1 / d1; y++ )
    {
      int C = -Math::intDiv ( e0 + b1 * y, a1 );
      for ( int x = C; x < C + m_omega * alpha0 / a1; x++ )
        {
          Vector2D X = H * Vector2D ( x, y );
          Vector2D I = calculate ( X );
	  Vector2D rem = calculateRemainder( X );

	  std::cout << "Rem = "<<rem[0]<<","<<rem[1]<<" "<<std::endl;

          setPavingRemainder ( I, rem, X );
        }
    }
}

/**
 * determines the set of pavings of the initial period - version using remains
 *
 */
void QAT::determinePavingsNoMultiply()
{
  Paving paving;

  int A = -Math::intDiv ( f0, d1 );
  int C = -Math::intDiv ( e0 + b1 * A, a1 );
  int Cp = Math::mod ( e0 + b1 * A, a1 );
  Vector2D X0 = H * Vector2D ( C, A );
  Vector2D I0 = calculate ( X0 );
  Vector2D I0p = m_matrix * X0 + m_vector - I0 * m_omega;

  int a10p = Math::intDiv ( a1, m_omega );
  int a11p = Math::mod ( a1, m_omega );
  int d10p = Math::intDiv ( d1, m_omega );
  int d11p = Math::mod ( d1, m_omega );

  int p = b1 - b10 * a1;
  int p0 = Math::intDiv ( p, m_omega );
  int p1 = Math::mod ( p, m_omega );

  for ( int y = A; y < A + m_omega * alpha1 / d1; y++ )
    {
      Vector2D X = X0;
      Vector2D I = I0;
      Vector2D Ip = I0p;
      for ( int x = C; x < C + m_omega * alpha0 / a1; x++ )
        {
          setPaving ( I, X );
          X += H[0];
          I[0] += a10p;
          Ip[0] += a11p;
          if ( Ip[0] >= m_omega )
            {
              I[0]++;
              Ip[0] -= m_omega;
            }
        }
      C -= b10;
      Cp += b11;
      X0 += H[1] - H0[0];
      I0[1] += d10p;
      I0p[1] += d11p;
      I0[0] += p0;
      I0p[0] += p1;
      if ( I0p[1] >= m_omega )
        {
          I0[1]++;
          I0p[1] -= m_omega;
        }
      if ( I0p[0] >= m_omega )
        {
          I0[0]++;
          I0p[0] -= m_omega;
        }
      if ( Cp >= a1 )
        {
          C--;
          Cp -= a1;
          X0 -= H[0];
          I0[0] -= a10p;
          I0p[0] -= a11p;
          if ( I0p[0] < 0 )
            {
              I0[0]--;
              I0p[0] += m_omega;
            }
        }
    }
}

/**
 * determines the set of pavings of the initial period - naive version
 *
 */
void QAT::determinePavingsByBoundingRect()
{
  for ( int i = 0; i < alpha0; i++ )
    for ( int j = 0; j < alpha1; j++ )
      setPaving ( i, j, determinePavingByBoundingRect ( i, j ) );
}

/**
 * computes the color of pixel P with floating point linear backward mapping method
 *
 * @param P the point to get color
 * @param image the source image
 *
 * @return the color
 */
const Color QAT::backwardColorLinear ( const Vector2D Pp, Image image ) const
  {
    float x = ( float ) Pp[0] / m_omega;
    float y = ( float ) Pp[1] / m_omega;
    float r = x - floor ( x );
    float l = 1 - r;
    float u = y - floor ( y );
    float d = 1 - u;
    Color color ( 0, 0, 0 );
    float count = 0;
    if ( 0 <= floor ( x ) && floor ( x ) < image.width() && 0 <= floor ( y ) && floor ( y ) < image.height() )
      {
        color += image.getColor ( floor ( x ), floor ( y ) ) * l * d;
        count += l * d;
      }
    if ( 0 <= floor ( x ) && floor ( x ) < image.width() && 0 <= ceil ( y ) && ceil ( y ) < image.height() )
      {
        color += image.getColor ( floor ( x ), ceil ( y ) ) * l * u;
        count += l * u;
      }
    if ( 0 <= ceil ( x ) && ceil ( x ) < image.width() && 0 <= floor ( y ) && floor ( y ) < image.height() )
      {
        color += image.getColor ( ceil ( x ), floor ( y ) ) * r * d;
        count += r * d;
      }
    if ( 0 <= ceil ( x ) && ceil ( x ) < image.width() && 0 <= ceil ( y ) && ceil ( y ) < image.height() )
      {
        color += image.getColor ( ceil ( x ), ceil ( y ) ) * r * u;
        count += r * u;
      }
    if ( count )
      color /= count;
    return color;
  }

/**
 * computes the color of pixel P with floating point Nearest Neighbor backward mapping method
 *
 * @param P the point to get color
 * @param image the source image
 *
 * @return the color
 */
const Color QAT::backwardColorNN ( const Vector2D Pp, Image image ) const
  {
    float x = ( float ) Pp[0] / m_omega;
    float y = ( float ) Pp[1] / m_omega;
    int rx = ( int ) rint ( x );
    int ry = ( int ) rint ( y );

    if ( rx < 0 )
      rx = 0;
    if ( rx >= image.width() )
      rx = image.width() - 1;
    if ( ry < 0 )
      ry = 0;
    if ( ry >= image.height() )
      ry = image.height() - 1;

    return image.getColor ( rx,ry );
  }

/**
 * applies the qat to an image
 *
 * @param image the initial image
 * @param backwardMapping use floating point backward mapping method (not compatible with any other option)
 * @param useBoundingRect true if using the naive method
 * @param usePeriodicity true if using the paving periodicity
 * @param noMultiply don't use many multiplications (not compatible with useBoundingRect)
 *
 * @return the final image
 */
const Image QAT::applyToImage ( Image image, InterpolationType interptype, bool useBoundingRect, bool usePeriodicity, bool noMultiply, bool fakeColor, bool is_inverse )
{
  int min_i, max_i, min_j, max_j;
  Image finalImage;
  if ( ( interptype == LINEAR ) or ( interptype == NN ) )
    {
      Matrix2x2 m_matrix_save;
      int m_omega_save;
      Vector2D m_vector_save;
      if(is_inverse)
	{
	  m_matrix_save = m_matrix;
	  m_omega_save = m_omega;
	  m_vector_save = m_vector;
	  inverse();
	}
      Matrix2x2 m = getImageBound ( image );
      min_i = m[0][0];
      min_j = m[1][0];
      max_i = m[0][1];
      max_j = m[1][1];
      finalImage = Image ( max_i - min_i + 1, max_j - min_j + 1 );
      if(is_inverse)
	{
	  m_matrix = m_matrix_save;
	  m_omega = m_omega_save;
	  m_vector = m_vector_save;
	}
      else
	inverse();
      Vector2D Pp = m_matrix * Vector2D ( min_i, min_j ) + m_vector;
      Vector2D IncrX = m_matrix * Vector2D ( 1, - ( max_j - min_j + 1 ) );
      Vector2D IncrY = m_matrix * Vector2D ( 0, 1 );
      if ( interptype == LINEAR )
        for ( int i = min_i; i <= max_i; i++ )
          {
            for ( int j = min_j; j <= max_j; j++ )
              {
                finalImage.setColor ( i - min_i, j - min_j, backwardColorLinear ( Pp, image ) );
                Pp += IncrY;
              }
            Pp += IncrX;
          }
      else
        for ( int i = min_i; i <= max_i; i++ )
          {
            for ( int j = min_j; j <= max_j; j++ )
              {
                finalImage.setColor ( i - min_i, j - min_j, backwardColorNN ( Pp, image ) );
                Pp += IncrY;
              }
            Pp += IncrX;
          }
    }
  else
    {
      Paving P,rP;
      int hShift, vShift;
      Matrix2x2 m_matrix_save;
      int m_omega_save;
      Vector2D m_vector_save;
      if(is_inverse)
	{
	  m_matrix_save = m_matrix;
	  m_omega_save = m_omega;
	  m_vector_save = m_vector;
	  inverse();
	}
      bool contracting = isContracting();

      if ( contracting )
        {
          std::cout << "Application contractante\n";
          Matrix2x2 m = getImageBound ( image );
          min_i = m[0][0];
          min_j = m[1][0];
          max_i = m[0][1];
          max_j = m[1][1];
          finalImage = Image ( max_i - min_i + 1, max_j - min_j + 1 );
        }
      else
        {
          int x, y;
          Matrix2x2 m = getImageBound ( image );
          hShift = m[0][0];
          vShift = m[1][0];
          x = m[0][1];
          y = m[1][1];
          if ( isInversible() )
            std::cout << "Inversible\n";
	  if(is_inverse)
	    {
	      m_matrix = m_matrix_save;
	      m_omega = m_omega_save;
	      m_vector = m_vector_save;
	    }
	  else
	    inverse();
          min_i = 0;
          max_i = image.width() - 1;
          min_j = 0;
          max_j = image.height() - 1;
          finalImage = Image ( x - hShift + 1, y - vShift + 1 );
        }
      compute();
      pavings = std::vector<Paving> ( alpha1 * alpha0, Paving() );
      pavingsRemainder = std::vector<Paving> ( alpha1 * alpha0, Paving() );
      Vector2D v ( 0, 0 );
      if ( usePeriodicity )
        if ( !useBoundingRect )
          if ( noMultiply )
            determinePavingsNoMultiply();
          else
            {
	      determinePavingsWithRemainders();
	    }
        else
          determinePavingsByBoundingRect();
      int j00 = Math::mod ( min_j, alpha1 );
      int v10 = Math::intDiv ( min_j, alpha1 );
      int i00 = Math::mod ( min_i + v10 * beta1, alpha0 );
      int v00 = Math::intDiv ( min_i + v10 * beta1, alpha0 );
      int i0, v0, j0, v1;
      int beta10 = Math::intDiv ( beta1, alpha0 );
      int beta11 = Math::mod ( beta1, alpha0 );
      for ( int i = min_i; i <= max_i; i++ )
        {
          if ( usePeriodicity && noMultiply )
            {
              i0 = i00;
              v0 = v00;
              j0 = j00;
              v1 = v10;
            }
          for ( int j = min_j; j <= max_j; j++ )
            {
              if ( usePeriodicity )
                {
                  if ( !noMultiply )
                    {
                      j0 = Math::mod ( j, alpha1 );
                      v1 = Math::intDiv ( j, alpha1 );
                      i0 = Math::mod ( i + v1 * beta1, alpha0 );
                      v0 = Math::intDiv ( i + v1 * beta1, alpha0 );
                    }
                  v = U0 * v0 + U1 * v1;
                  P = getPaving ( i0, j0 );
		  rP = getPavingRemainder ( i0, j0 );
                }
              else
                if ( useBoundingRect )
                  P = determinePavingByBoundingRect ( i, j );
                else
                  if ( noMultiply )
                    P = determinePavingNoMultiply ( i, j );
                  else
                    {
		      P = determinePaving ( i, j );
		    }
              if ( contracting )
                finalImage.setColor ( i - min_i, j - min_j, image.colorOfPaving ( P, v ) );
              else
                if ( fakeColor )
                  finalImage.colorizePavingRemainder ( P, rP, Vector2D ( -hShift, -vShift ) + v, 
					      image.getFakeColor ( i, j ) ,
					      image.getFakeColor ( i, j ), 
					      image.getFakeColor ( i, j ),
					      image.getFakeColor ( i, j ),
					      image.getFakeColor ( i, j ),m_omega);
                else
                  finalImage.colorizePavingRemainder ( P, rP, Vector2D ( -hShift, -vShift ) + v, 
					      image.getColor ( i, j ), 
					      image.getColor ( i, j ), 
					      image.getColor ( i, j ),
					      image.getColor ( i, j ),
					      image.getColor ( i, j ),m_omega);

              if ( usePeriodicity && noMultiply )
                {
                  j0++;
                  if ( j0 == alpha1 )
                    {
                      v1++;
                      j0 = 0;
                      v0 += beta10;
                      i0 += beta11;
                      if ( i0 >= alpha0 )
                        {
                          v0++;
                          i0 -= alpha0;
                        }
                    }
                }
            }
          if ( usePeriodicity && noMultiply )
            {
              i00++;
              if ( i00 == alpha0 )
                {
                  v00++;
                  i00 = 0;
                }
            }
        }
    }

  return finalImage ;
}

/**
 * inverses the qat
 *
 */
void QAT::inverse()
{

  //m_matrix.print();
  //Matrix2x2 tmp = m_matrix;
  int det = m_matrix.det();
  a0 = m_matrix[1][1];
  b0 = -m_matrix[1][0];
  c0 = -m_matrix[0][1];
  d0 = m_matrix[0][0];
  m_matrix[0][0] = a0 * m_omega;
  m_matrix[1][0] = b0 * m_omega;
  m_matrix[0][1] = c0 * m_omega;
  m_matrix[1][1] = d0 * m_omega;
  if ( det < 0 )
    m_matrix *= -1;
  m_vector = -m_matrix * m_vector;
  m_omega = abs ( det );
// m_matrix.print();
//tmp = tmp*m_matrix;
//std::cout <<"Id =";tmp.print();

}

/**
 * computes the minimum and maximum coordinates in the final image : min_x, min_y, max_x, max_y
 *
 * @param image the initial image
 *
 * @return the matrix : (min_i, min_j, max_i, max_j)
 */
const Matrix2x2 QAT::getImageBound ( const Image image )
{
  Vector2D S1 = calculate ( Vector2D ( 0, 0 ) );
  Vector2D S2 = calculate ( Vector2D ( 0, image.height() ) );
  Vector2D S3 = calculate ( Vector2D ( image.width(), image.height() ) );
  Vector2D S4 = calculate ( Vector2D ( image.width(), 0 ) );

  return Matrix2x2 ( std::min ( std::min ( S1[0], S2[0] ), std::min ( S3[0], S4[0] ) ),
                     std::min ( std::min ( S1[1], S2[1] ), std::min ( S3[1], S4[1] ) ),
                     std::max ( std::max ( S1[0], S2[0] ), std::max ( S3[0], S4[0] ) ),
                     std::max ( std::max ( S1[1], S2[1] ), std::max ( S3[1], S4[1] ) ) );

}

/**
 * computes the necessary constants
 *
 */
void QAT::compute()
{
  int cp0, dp0, u0, v0;
  Math::couple v;

  a0 = m_matrix[0][0];
  b0 = m_matrix[1][0];
  c0 = m_matrix[0][1];
  d0 = m_matrix[1][1];
  e0 = m_vector[0];
  f0 = m_vector[1];

  if ( c0 != 0 )
    {
      v = Math::extendedEuclidean ( c0, d0 );
      cp0 = c0 / Math::pgcd ( c0, d0 );
      dp0 = d0 / Math::pgcd ( c0, d0 );
      u0 = v[0];
      v0 = v[1];

      a1 = a0 * dp0 - b0 * cp0;
      b1 = a0 * u0 + b0 * v0;
      c1 = 0;
      d1 = c0 * u0 + d0 * v0;
    }
  else
    {
      cp0 = 0;
      dp0 = 1;
      u0 = 0;
      v0 = 1;

      a1 = a0;
      b1 = b0;
      c1 = c0;
      d1 = d0;
    }

  H = Matrix2x2 ( dp0, u0, -cp0, v0 );

  if ( a1 < 0 )
    {
      a1 = -a1;
      H *= Matrix2x2 ( -1, 0, 0, 1 );
    }

  b10 = Math::intDiv ( b1, a1 );
  b11 = Math::mod ( b1, a1 );

  H0[0] = H[0] * b10;

  alpha0 = a1 / Math::pgcd ( a1, m_omega );
  U0 = H * Vector2D ( m_omega / Math::pgcd ( a1, m_omega ), 0 );

  {
    int d1p = d1 / Math::pgcd ( d1, m_omega );
    int omegap = m_omega / Math::pgcd ( d1, m_omega );
    int a1p = a1 / Math::pgcd ( Math::pgcd ( a1, m_omega ), b1 * omegap );
    int phi = b1 * omegap / Math::pgcd ( Math::pgcd ( a1, m_omega ), b1 * omegap );
    int omegas = m_omega / Math::pgcd ( Math::pgcd ( a1, m_omega ), b1 * omegap );
    Math::couple vect = Math::extendedEuclidean ( a1p, omegas );
    int u1 = vect[0];
    int v1 = vect[1];
    alpha1 = d1p * Math::pgcd ( a1p, omegas );
    beta1 = -phi * v1;
    U1 = H * Vector2D ( -phi * u1, omegap * Math::pgcd ( a1p, omegas ) );
  }
}

/**
 * computes the image of a point
 *
 * @param p the initial point
 *
 * @return the image of the point
 */
const Vector2D QAT::calculate ( const Vector2D p ) const
  {
    return ( m_matrix * p + m_vector ) / m_omega;
  }

/**
* computes the remainder of a point
*
* @param p the initial point
*
* @return the image of the point
*/
const Vector2D QAT::calculateRemainder ( const Vector2D p ) const
{
  return ( m_matrix * p + m_vector ) % m_omega;
}


