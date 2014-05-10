/**
 * @file   qat.cpp
 * @author Valentin Blot <valentin@Valentin-E21ter>
 * @date   Thu Jun 19 14:28:51 2008
 *
 * @brief  quasi-affine transforms utilities
 *
 *
 */

#include "qat3d.hpp"

/**
 * default constructor
 *
 */
QAT3D::QAT3D()
{
}

/**
 * construct from a matrix, an integer, and a vector
 *
 * @param matrix the matrix of the qat
 * @param omega the omega coefficient
 * @param vect the translation vector of the qat
 */
QAT3D::QAT3D ( Matrix3x3 matrix, int omega, Vector3D vect )
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
bool QAT3D::isContracting() const
{
  return m_omega * m_omega * m_omega > abs ( m_matrix.det() );
}

/**
 * determines the set of antecedents of a point by the qat - version with multiplications
 *
 * @param i X-axis coordinate of the point
 * @param j Y-axis coordinate of the point
 * @param k Z-axis coordinate of the point
 *
 * @return the paving containing the antecedents
 */
const Paving3D QAT3D::determinePaving ( int i, int j, int k ) const
{
  Paving3D paving;
  
  int A = -Math::intDiv ( -k * m_omega + l0, i1 );
  int B = -Math::intDiv ( - ( k + 1 ) * m_omega + l0, i1 );
  
  for ( int z = A; z < B; z++ )
  {
    int C = -Math::intDiv ( -j * m_omega + k0 + f1 * z, e1 );
    int D = -Math::intDiv ( - ( j + 1 ) * m_omega + k0 + f1 * z, e1 );
    for ( int y = C; y < D; y++ )
    {
      int E = -Math::intDiv ( -i * m_omega + j0 + c1 * z + b1 * y, a1 );
      int F = -Math::intDiv ( - ( i + 1 ) * m_omega + j0 + c1 * z + b1 * y, a1 );
      for ( int x = E; x < F; x++ )
        paving.addPoint ( H * Vector3D ( x, y, z ) );
    }
  }
  return paving;
}

/**
 * determines the set of antecedents of a point by the qat - version using remains
 *
 * @param i X-axis coordinate of the point
 * @param j Y-axis coordinate of the point
 * @param k Z-axis coordinate of the point
 *
 * @return the paving containing the antecedents
 */
const Paving3D QAT3D::determinePavingNoMultiply ( int i, int j, int k ) const
{
  Paving3D paving;
  
  int A = -Math::intDiv ( -k * m_omega + l0, i1 );
  int B = -Math::intDiv ( - ( k + 1 ) * m_omega + l0, i1 );
  int C = -Math::intDiv ( -j * m_omega + k0 + f1 * A, e1 );
  int D = -Math::intDiv ( - ( j + 1 ) * m_omega + k0 + f1 * A, e1 );
  int E0 = -Math::intDiv ( -i * m_omega + j0 + c1 * A + b1 * C, a1 );
  int F0 = -Math::intDiv ( - ( i + 1 ) * m_omega + j0 + c1 * A + b1 * C, a1 );
  int Cp = -Math::mod ( -j * m_omega + k0 + f1 * A, e1 );
  int Dp = -Math::mod ( - ( j + 1 ) * m_omega + k0 + f1 * A, e1 );
  int E0p = -Math::mod ( -i * m_omega + j0 + c1 * A + b1 * C, a1 );
  int F0p = -Math::mod ( - ( i + 1 ) * m_omega + j0 + c1 * A + b1 * C, a1 );
  Vector3D X0 = H * Vector3D ( E0, C, A );
  
  for ( int z = A; z < B; z++ )
  {
    int E = E0, F = F0, Ep = E0p, Fp = F0p;
    Vector3D X1 = X0;
    for ( int y = C; y < D; y++ )
    {
      Vector3D X2 = X1;
      for ( int x = E; x < F; x++ )
      {
        paving.addPoint ( X2 );
        X2 += H[0];
      }
      E -= b10;
      Ep += b11;
      X1 += H[1] - H0[0];
      if ( Ep >= a1 )
      {
        E--;
        Ep -= a1;
        X1 -= H[0];
      }
      F -= b10;
      Fp += b11;
      if ( Fp >= a1 )
      {
        F--;
        Fp -= a1;
      }
    }
    C -= f10;
    Cp += f11;
    X0 += H[2] - H0[1];
    if ( Cp >= e1 )
    {
      C--;
      Cp -= e1;
      E0 -= p0;
      E0p += p1;
      X0 -= H[1] + H1[0];
      F0 -= p0;
      F0p += p1;
    }
    else
    {
      E0 -= q0;
      E0p += q1;
      X0 -= H2[0];
      F0 -= q0;
      F0p += q1;
    }
    if ( E0p >= a1 )
    {
      E0--;
      E0p -= a1;
      X0 -= H[0];
    }
    if ( F0p >= a1 )
    {
      F0--;
      F0p -= a1;
    }
    D -= f10;
    Dp += f11;
    if ( Dp >= e1 )
    {
      D--;
      Dp -= e1;
    }
  }
  return paving;
}

/**
 * determines the set of antecedents of a point by the qat - naive version
 *
 * @param i X-axis coordinate of the point
 * @param j Y-axis coordinate of the point
 * @param k Z-axis coordinate of the point
 *
 * @return the paving containing the antecedents
 */
const Paving3D QAT3D::determinePavingByBoundingRect ( int i, int j, int k ) const
{
  Paving3D paving;
  int delta = m_matrix.det();
  for ( int U = m_omega * i; U < m_omega * ( i + 1 ); U++ )
    for ( int V = m_omega * j; V < m_omega * ( j + 1 ); V++ )
      for ( int W = m_omega * k; W < m_omega * ( k + 1 ); W++ )
      {
    int x = ( e0 * i0 - h0 * f0 ) * U + ( h0 * c0 - b0 * i0 ) * V + ( b0 * f0 - e0 * c0 ) * W;
    int y = ( g0 * f0 - d0 * i0 ) * U + ( a0 * i0 - g0 * c0 ) * V + ( d0 * c0 - a0 * f0 ) * W;
    int z = ( d0 * h0 - g0 * e0 ) * U + ( g0 * b0 - a0 * h0 ) * V + ( a0 * e0 - d0 * b0 ) * W;
    if ( x % delta == 0 && y % delta == 0 && z % delta == 0 )
    {
      x /= delta;
      y /= delta;
      z /= delta;
      paving.addPoint ( Vector3D ( x, y, z ) );
    }
  }
  
  return paving;
}

/**
 * adds a point to a paving in the set of pavings of the initial period
 *
 * @param I coordinates of the paving
 * @param P the point to add
 */
void QAT3D::setPaving ( const Vector3D I, const Vector3D P )
{
  if ( 0 <= I[0] && I[0] < alpha0 && 0 <= I[1] && I[1] < alpha1 && 0 <= I[2] && I[2] < alpha2 )
    pavings[I[2] * alpha1 * alpha0 + I[1] * alpha0 + I[0]].addPoint ( P );
}

/**
 * sets a paving in the set of pavings of the initial period
 *
 * @param i first coordinate of the paving
 * @param j second coordinate of the paving
 * @param k third coordinate of the paving
 * @param P the paving
 */
void QAT3D::setPaving ( int i, int j, int k, const Paving3D P )
{
  if ( 0 <= i && i < alpha0 && 0 <= j && j < alpha1 && 0 <= k && k < alpha2 )
    pavings[k * alpha1 * alpha0 + j * alpha0 + i] = P;
}

/**
 * gets a paving in the set of pavings of the initial period
 *
 * @param i first coordinate of the paving
 * @param j second coordinate of the paving
 * @param k third coordinate of the paving
 *
 * @return the paving
 */
const Paving3D QAT3D::getPaving ( int i, int j, int k ) const
{
  if ( 0 <= i && i < alpha0 && 0 <= j && j < alpha1 && 0 <= k && k < alpha2 )
    return pavings[k * alpha1 * alpha0 + j * alpha0 + i];
  return Paving3D();
}

/**
 * determines the set of pavings of the initial period - version with multiplications
 *
 */
void QAT3D::determinePavings()
{
  int A = -Math::intDiv ( l0, i1 );
  for ( int z = A; z < A + m_omega * alpha2 / i1; z++ )
  {
    int C = -Math::intDiv ( k0 + f1 * z, e1 );
    for ( int y = C; y < C + m_omega * alpha1 / e1; y++ )
    {
      int E = -Math::intDiv ( j0 + c1 * z + b1 * y, a1 );
      for ( int x = E ; x < E + m_omega * alpha0 / a1; x++ )
      {
        Vector3D X = H * Vector3D ( x, y, z );
        Vector3D I = calculate ( X );
        setPaving ( I, X );
      }
    }
  }
}

/**
 * determines the set of pavings of the initial period - version using remains
 *
 */
void QAT3D::determinePavingsNoMultiply()
{
  // TODO
}

/**
 * determines the set of pavings of the initial period - naive version
 *
 */
void QAT3D::determinePavingsByBoundingRect()
{
  for ( int i = 0; i < alpha0; i++ )
    for ( int j = 0; j < alpha1; j++ )
      for ( int k = 0; k < alpha2; k++ )
        setPaving ( i, j, k, determinePavingByBoundingRect ( i, j, k ) );
}

/**
 * computes the color of pixel P with floating point linear backward mapping method
 *
 * @param P the point to get color
 * @param image the source image
 *
 * @return the color
 */
const Color QAT3D::backwardColorLinear ( const Vector3D Pp, Image3D image ) const
{
  double x = ( double ) Pp[0] / m_omega;
  double y = ( double ) Pp[1] / m_omega;
  double z = ( double ) Pp[2] / m_omega;
  double r = x - floor ( x );
  double l = 1 - r;
  double u = y - floor ( y );
  double d = 1 - u;
  double b = z - floor ( z );
  double f = 1 - b;
  Color color ( 0, 0, 0 );
  double count = 0;
  if ( 0 <= floor ( x ) && floor ( x ) < image.width() )
  {
    if ( 0 <= floor ( y ) && floor ( y ) < image.height() )
    {
      if ( 0 <= floor ( z ) && floor ( z ) < image.depth() )
      {
        color += image.getColor ( floor ( x ), floor ( y ), floor ( z ) ) * l * d * f;
        count += l * d * f;
      }
      if ( 0 <= ceil ( z ) && ceil ( z ) < image.depth() )
      {
        color += image.getColor ( floor ( x ), floor ( y ), ceil ( z ) ) * l * d * b;
        count += l * d * b;
      }
    }
    if ( 0 <= ceil ( y ) && ceil ( y ) < image.height() )
    {
      if ( 0 <= floor ( z ) && floor ( z ) < image.depth() )
      {
        color += image.getColor ( floor ( x ), ceil ( y ), floor ( z ) ) * l * u * f;
        count += l * u * f;
      }
      if ( 0 <= ceil ( z ) && ceil ( z ) < image.depth() )
      {
        color += image.getColor ( floor ( x ), ceil ( y ), ceil ( z ) ) * l * u * b;
        count += l * u * b;
      }
    }
  }
  if ( 0 <= ceil ( x ) && ceil ( x ) < image.width() )
  {
    if ( 0 <= floor ( y ) && floor ( y ) < image.height() )
    {
      if ( 0 <= floor ( z ) && floor ( z ) < image.depth() )
      {
        color += image.getColor ( ceil ( x ), floor ( y ), floor ( z ) ) * r * d * f;
        count += r * d * f;
      }
      if ( 0 <= ceil ( z ) && ceil ( z ) < image.depth() )
      {
        color += image.getColor ( ceil ( x ), floor ( y ), ceil ( z ) ) * r * d * b;
        count += r * d * b;
      }
    }
    if ( 0 <= ceil ( y ) && ceil ( y ) < image.height() )
    {
      if ( 0 <= floor ( z ) && floor ( z ) < image.depth() )
      {
        color += image.getColor ( ceil ( x ), ceil ( y ), floor ( z ) ) * r * u * f;
        count += r * u * f;
      }
      if ( 0 <= ceil ( z ) && ceil ( z ) < image.depth() )
      {
        color += image.getColor ( ceil ( x ), ceil ( y ), ceil ( z ) ) * r * u * b;
        count += r * u * b;
      }
    }
  }
  if ( count )
    color /= count;
  return color;
}

/**
 * computes the color of pixel P with floating point nearest neighbor backward mapping method
 *
 * @param P the point to get color
 * @param image the source image
 *
 * @return the color
 */
const Color QAT3D::backwardColorNN ( const Vector3D Pp, Image3D image ) const
{
  int rx = ( int ) rint ( ( double ) Pp[0] / m_omega );
  int ry = ( int ) rint ( ( double ) Pp[1] / m_omega );
  int rz = ( int ) rint ( ( double ) Pp[2] / m_omega );
  
  if ( rz < 0 ) rz = 0;
  if ( rx < 0 ) rx = 0;
  if ( ry < 0 ) ry = 0;
  if ( rx >= image.width() ) rx = image.width();
  if ( ry >= image.height() ) ry = image.height();
  if ( rz >= image.depth() ) rz = image.depth();
  
  return image.getColor ( rx,ry,rz );
  
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
const Image3D QAT3D::applyToImage ( const Image3D image, InterpolationType interp, bool useBoundingRect, bool usePeriodicity, bool noMultiply, bool fakeColor, bool is_inverse )
{
  //is_inverse = false;
  int min_i, max_i, min_j, max_j, min_k, max_k;
  Image3D finalImage;
  if ( ( interp == LINEAR ) or ( interp == NN ) )
  {
    Matrix3x3 m_matrix_save;
    int m_omega_save;
    Vector3D m_vector_save;
    if(is_inverse)
    {
      m_matrix_save = m_matrix;
      m_omega_save = m_omega;
      m_vector_save = m_vector;
      inverse();
    }
    Matrix3x3 m = getImageBound ( image );
    min_i = m[0][0];
    min_j = m[1][0];
    min_k = m[2][0];
    max_i = m[0][1];
    max_j = m[1][1];
    max_k = m[2][1];
    finalImage = Image3D ( max_i - min_i + 1, max_j - min_j + 1, max_k - min_k + 1 );
    std::cout << "Transformed image size: "<<max_i-min_i+1<<"x"<<max_j - min_j + 1<<"x"<<max_k - min_k + 1 <<std::endl;
    if(is_inverse)
    {
      m_matrix = m_matrix_save;
      m_omega = m_omega_save;
      m_vector = m_vector_save;
    }
    else
      inverse();
    Vector3D Pp = m_matrix * Vector3D ( min_i, min_j, min_k ) + m_vector;
    Vector3D IncrX = m_matrix * Vector3D ( 1, - ( max_j - min_j + 1 ), 0 );
    Vector3D IncrY = m_matrix * Vector3D ( 0, 1, - ( max_k - min_k + 1 ) );
    Vector3D IncrZ = m_matrix * Vector3D ( 0, 0, 1 );
    if ( interp == LINEAR )
      for ( int i = min_i; i <= max_i; i++ )
      {
      std::cout << "i="<<i-min_i<<"/"<<max_i-min_i<<std::endl;
      for ( int j = min_j; j <= max_j; j++ )
      {
        std::cout << "j="<<j-min_j<<"/"<<max_j-min_j<<std::endl;
        for ( int k = min_k; k <= max_k; k++ )
        {
          std::cout << "k="<<k-min_k<<"/"<<max_k-min_k<<std::endl;
          finalImage.setColor ( i - min_i, j - min_j, k - min_k, backwardColorLinear ( Pp, image ) );
          Pp += IncrZ;
        }
        Pp += IncrY;
      }
      Pp += IncrX;
    }
    else
      for ( int i = min_i; i <= max_i; i++ )
      {
      for ( int j = min_j; j <= max_j; j++ )
      {
        for ( int k = min_k; k <= max_k; k++ )
        {
          finalImage.setColor ( i - min_i, j - min_j, k - min_k, backwardColorNN ( Pp, image ) );
          Pp += IncrZ;
        }
        Pp += IncrY;
      }
      Pp += IncrX;
    }
  }
  else
  {
    QAT3D h;
    Paving3D P;
    int hShift, vShift, pShift;
    Matrix3x3 m_matrix_save;
    int m_omega_save;
    Vector3D m_vector_save;
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
      Matrix3x3 m = getImageBound ( image );
      min_i = m[0][0];
      min_j = m[1][0];
      min_k = m[2][0];
      max_i = m[0][1];
      max_j = m[1][1];
      max_k = m[2][1];
      finalImage = Image3D ( max_i - min_i + 1, max_j - min_j + 1, max_k - min_k + 1 );
      std::cout << "Transformed image size: "<<max_i-min_i+1<<"x"<<max_j - min_j + 1<<"x"<<max_k - min_k + 1 <<std::endl;
    }
    else
    {
      int x, y, z;
      Matrix3x3 m = getImageBound ( image );
      hShift = m[0][0];
      vShift = m[1][0];
      pShift = m[2][0];
      x = m[0][1];
      y = m[1][1];
      z = m[2][1];
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
      min_k = 0;
      max_k = image.depth() - 1;
      finalImage = Image3D ( x - hShift + 1, y - vShift + 1, z - pShift + 1 );
      std::cout << "Transformed image size: "<<x - hShift + 1<<"x"<<y - vShift + 1<<"x"<<z - pShift + 1<<std::endl;
    }
    compute();
    pavings = std::vector<Paving3D> ( alpha0 * alpha1 * alpha2, Paving3D() );
    Vector3D v ( 0, 0, 0 );
    if ( usePeriodicity )
    {
      if ( !useBoundingRect )
      {
        if ( noMultiply )
          determinePavingsNoMultiply();
        else
          determinePavings();
      }
      else
        determinePavingsByBoundingRect();
    }

    for ( int i = min_i; i <= max_i; i++ )
    {
      std::cout << ( i - min_i ) / ( max_i - min_i ) << "%\r";
      if ( usePeriodicity && noMultiply )
      {
        std::cerr << "Noo_multiply not yet implemented in 3D"<<std::endl;
        // TODO
      }
      for ( int j = min_j; j <= max_j; j++ )
      {
        if ( usePeriodicity && noMultiply )
        {
          std::cerr << "Noo_multiply not yet implemented in 3D"<<std::endl;
          // TODO
        }
        for ( int k = min_k; k <= max_k; k++ )
        {
          if ( usePeriodicity )
          {
            if ( !noMultiply )
            {
              int kp = Math::mod ( k, alpha2 );
              int ks = Math::intDiv ( k, alpha2 );
              int jp = Math::mod ( j + ks * gamma2, alpha1 );
              int js = Math::intDiv ( j + ks * gamma2, alpha1 );
              int ip = Math::mod ( i + ks * beta2 + js * beta1, alpha0 );
              int is = Math::intDiv ( i + ks * beta2 + js * beta1, alpha0 );
              v = U0 * is + U1 * js + U2 * ks;
              P = getPaving ( ip, jp, kp );
            }
            else
            {
              std::cerr << "No_multiply not yet implemented in 3D"<<std::endl;
              // TODO
            }
          }
          else
          {
            if ( useBoundingRect )
            {
              P = determinePavingByBoundingRect ( i, j, k );
            }
            else
            {
              if ( noMultiply )
                P = determinePavingNoMultiply ( i, j, k );
              else
                P = determinePaving ( i, j, k );
            }
          }

          if ( contracting )
            finalImage.setColor ( i - min_i, j - min_j, k - min_k, image.colorOfPaving ( P, v ) );
          else
          {
            if ( fakeColor )
              finalImage.colorizePaving ( P, Vector3D ( -hShift, -vShift, -pShift ) + v, image.getFakeColor ( i, j, k ) );
            else
              finalImage.colorizePaving ( P, Vector3D ( -hShift, -vShift, -pShift ) + v, image.getColor ( i, j, k ) );
          }

          if ( usePeriodicity && noMultiply )
          {
            std::cerr << "Noo_multiply not yet implemented in 3D"<<std::endl;
            // TODO
          }
        }
        if ( usePeriodicity && noMultiply )
        {
          std::cerr << "Noo_multiply not yet implemented in 3D"<<std::endl;
          // TODO
        }
      }
      if ( usePeriodicity && noMultiply )
      {
        // TODO
        std::cerr << "Noo_multiply not yet implemented in 3D"<<std::endl;
      }
    }
  }
  
  return finalImage ;
}

/**
 * inverses the qat
 *
 */
void QAT3D::inverse()
{
  int det = m_matrix.det();
  a0 = m_matrix[1][1] * m_matrix[2][2] - m_matrix[2][1] * m_matrix[1][2];
  b0 = m_matrix[2][1] * m_matrix[0][2] - m_matrix[0][1] * m_matrix[2][2];
  c0 = m_matrix[0][1] * m_matrix[1][2] - m_matrix[1][1] * m_matrix[0][2];
  d0 = m_matrix[2][0] * m_matrix[1][2] - m_matrix[1][0] * m_matrix[2][2];
  e0 = m_matrix[0][0] * m_matrix[2][2] - m_matrix[2][0] * m_matrix[0][2];
  f0 = m_matrix[1][0] * m_matrix[0][2] - m_matrix[0][0] * m_matrix[1][2];
  g0 = m_matrix[1][0] * m_matrix[2][1] - m_matrix[2][0] * m_matrix[1][1];
  h0 = m_matrix[2][0] * m_matrix[0][1] - m_matrix[0][0] * m_matrix[2][1];
  i0 = m_matrix[0][0] * m_matrix[1][1] - m_matrix[1][0] * m_matrix[0][1];
  m_matrix[0][0] = a0 * m_omega;
  m_matrix[0][1] = b0 * m_omega;
  m_matrix[0][2] = c0 * m_omega;
  m_matrix[1][0] = d0 * m_omega;
  m_matrix[1][1] = e0 * m_omega;
  m_matrix[1][2] = f0 * m_omega;
  m_matrix[2][0] = g0 * m_omega;
  m_matrix[2][1] = h0 * m_omega;
  m_matrix[2][2] = i0 * m_omega;
  if ( det < 0 )
    m_matrix *= -1;
  m_vector = -m_matrix * m_vector;
  m_omega = abs ( det );
}

/**
 * computes the minimum and maximum coordinates in the final image : min_x, min_y, min_z, max_x, max_y, max_z
 *
 * @param image the initial image
 *
 * @return the matrix : (min_i, min_j, min_k, max_i, max_j, max_k, 0, 0, 0)
 */
const Matrix3x3 QAT3D::getImageBound ( const Image3D image )
{
  Vector3D S1 = calculate ( Vector3D ( 0, 0, 0 ) );
  Vector3D S2 = calculate ( Vector3D ( image.width() - 1, 0, 0 ) );
  Vector3D S3 = calculate ( Vector3D ( 0, image.height() - 1, 0 ) );
  Vector3D S4 = calculate ( Vector3D ( 0, 0, image.depth() - 1 ) );
  Vector3D S5 = calculate ( Vector3D ( image.width() - 1, image.height() - 1, 0 ) );
  Vector3D S6 = calculate ( Vector3D ( image.width() - 1, 0, image.depth() - 1 ) );
  Vector3D S7 = calculate ( Vector3D ( 0, image.height() - 1, image.depth() - 1 ) );
  Vector3D S8 = calculate ( Vector3D ( image.width() - 1, image.height() - 1, image.depth() - 1 ) );
  return Matrix3x3 ( std::min ( std::min ( std::min ( S1[0], S2[0] ), std::min ( S3[0], S4[0] ) ), std::min ( std::min ( S5[0], S6[0] ), std::min ( S7[0], S8[0] ) ) ),
                     std::min ( std::min ( std::min ( S1[1], S2[1] ), std::min ( S3[1], S4[1] ) ), std::min ( std::min ( S5[1], S6[1] ), std::min ( S7[1], S8[1] ) ) ),
                     std::min ( std::min ( std::min ( S1[2], S2[2] ), std::min ( S3[2], S4[2] ) ), std::min ( std::min ( S5[2], S6[2] ), std::min ( S7[2], S8[2] ) ) ),
                     std::max ( std::max ( std::max ( S1[0], S2[0] ), std::max ( S3[0], S4[0] ) ), std::max ( std::max ( S5[0], S6[0] ), std::max ( S7[0], S8[0] ) ) ),
                     std::max ( std::max ( std::max ( S1[1], S2[1] ), std::max ( S3[1], S4[1] ) ), std::max ( std::max ( S5[1], S6[1] ), std::max ( S7[1], S8[1] ) ) ),
                     std::max ( std::max ( std::max ( S1[2], S2[2] ), std::max ( S3[2], S4[2] ) ), std::max ( std::max ( S5[2], S6[2] ), std::max ( S7[2], S8[2] ) ) ),
                     0,
                     0,
                     0 );
}

/**
 * computes the necessary constants
 *
 */
void QAT3D::compute()
{
  Math::couple v;
  Matrix3x3 m;
  int a2, b2, c2, d2, e2, f2, g2, h2, i2;
  int a3, b3, c3, d3, e3, f3, g3, h3, i3;
  int gp0, hp0, u0, v0;
  int hp1, ip1, u1, v1;
  int dp2, ep2, u2, v2;
  
  a0 = m_matrix[0][0];
  b0 = m_matrix[1][0];
  c0 = m_matrix[2][0];
  d0 = m_matrix[0][1];
  e0 = m_matrix[1][1];
  f0 = m_matrix[2][1];
  g0 = m_matrix[0][2];
  h0 = m_matrix[1][2];
  i0 = m_matrix[2][2];
  
  j0 = m_vector[0];
  k0 = m_vector[1];
  l0 = m_vector[2];
  
  if ( g0 != 0 )
  {
    v = Math::extendedEuclidean ( g0, h0 );
    gp0 = g0 / Math::pgcd ( g0, h0 );
    hp0 = h0 / Math::pgcd ( g0, h0 );
    u0 = v[0];
    v0 = v[1];
    
    a1 = hp0 * a0 - gp0 * b0;
    b1 = u0 * a0 + v0 * b0;
    c1 = c0;
    d1 = hp0 * d0 - gp0 * e0;
    e1 = u0 * d0 + v0 * e0;
    f1 = f0;
    g1 = 0;
    h1 = u0 * g0 + v0 * h0;
    i1 = i0;
  }
  else
  {
    gp0 = 0;
    hp0 = 1;
    u0 = 0;
    v0 = 1;
    
    a1 = a0;
    b1 = b0;
    c1 = c0;
    d1 = d0;
    e1 = e0;
    f1 = f0;
    g1 = g0;
    h1 = h0;
    i1 = i0;
  }
  
  if ( h1 != 0 )
  {
    v = Math::extendedEuclidean ( h1, i1 );
    hp1 = h1 / Math::pgcd ( h1, i1 );
    ip1 = i1 / Math::pgcd ( h1, i1 );
    u1 = v[0];
    v1 = v[1];
    
    a2 = a1;
    b2 = ip1 * b1 - hp1 * c1;
    c2 = u1 * b1 + v1 * c1;
    d2 = d1;
    e2 = ip1 * e1 - hp1 * f1;
    f2 = u1 * e1 + v1 * f1;
    g2 = 0;
    h2 = 0;
    i2 = u1 * h1 + v1 * i1;
  }
  else
  {
    hp1 = 0;
    ip1 = 1;
    u1 = 0;
    v1 = 1;
    
    a2 = a1;
    b2 = b1;
    c2 = c1;
    d2 = d1;
    e2 = e1;
    f2 = f1;
    g2 = g1;
    h2 = h1;
    i2 = i1;
  }
  
  if ( d2 != 0 )
  {
    v = Math::extendedEuclidean ( d2, e2 );
    dp2 = d2 / Math::pgcd ( d2, e2 );
    ep2 = e2 / Math::pgcd ( d2, e2 );
    u2 = v[0];
    v2 = v[1];
    
    a3 = ep2 * a2 - dp2 * b2;
    b3 = u2 * a2 + v2 * b2;
    c3 = c2;
    d3 = 0;
    e3 = u2 * d2 + v2 * e2;
    f3 = f2;
    g3 = 0;
    h3 = 0;
    i3 = i2;
  }
  else
  {
    dp2 = 0;
    ep2 = 1;
    u2 = 0;
    v2 = 1;
    
    a3 = a2;
    b3 = b2;
    c3 = c2;
    d3 = d2;
    e3 = e2;
    f3 = f2;
    g3 = g2;
    h3 = h2;
    i3 = i2;
  }
  
  H = Matrix3x3 ( ep2 * hp0 - dp2 * ip1 * u0, u2 * hp0 + v2 * ip1 * u0, u0 * u1,
                  -ep2 * gp0 - dp2 * ip1 * v0, -u2 * gp0 + v2 * ip1 * v0, v0 * u1,
                  dp2 * hp1, -hp1 * v2, v1 );
  
  a1 = a3;
  b1 = b3;
  c1 = c3;
  d1 = d3;
  e1 = e3;
  f1 = f3;
  g1 = g3;
  h1 = h3;
  i1 = i3;
  
  if ( a1 < 0 )
  {
    a1 = -a1;
    H *= Matrix3x3 ( -1, 0, 0, 0, 1, 0, 0, 0, 1 );
  }
  if ( e1 < 0 )
  {
    b1 = -b1;
    e1 = -e1;
    H *= Matrix3x3 ( 1, 0, 0, 0, -1, 0, 0, 0, 1 );
  }
  
  f10 = Math::intDiv ( f1, e1 );
  f11 = Math::mod ( f1, e1 );
  
  b10 = Math::intDiv ( b1, a1 );
  b11 = Math::mod ( b1, a1 );
  
  p = -b1 * ( f10 + 1 ) + c1;
  p0 = Math::intDiv ( p, a1 );
  p1 = Math::mod ( p, a1 );
  
  q = -b1 * f10 + c1;
  q0 = Math::intDiv ( q, a1 );
  q1 = Math::mod ( q, a1 );
  
  H0[0] = H[0] * b10;
  H0[1] = H[1] * f10;
  H1[0] = H[0] * p0;
  H2[0] = H[0] * q0;
  
  alpha0 = a1 / Math::pgcd ( a1, m_omega );
  U0 = H * Vector3D ( m_omega / Math::pgcd ( a1, m_omega ), 0, 0 );
  
  {
    int e1p = e1 / Math::pgcd ( e1, m_omega );
    int omegap = m_omega / Math::pgcd ( e1, m_omega );
    int a1p = a1 / Math::pgcd ( Math::pgcd ( m_omega, a1 ), b1 * omegap );
    int phi = b1 * omegap / Math::pgcd ( Math::pgcd ( m_omega, a1 ), b1 * omegap );
    int omegas = m_omega / Math::pgcd ( Math::pgcd ( m_omega, a1 ), b1 * omegap );
    Math::couple vect = Math::extendedEuclidean ( a1p, omegas );
    int u1 = vect[0];
    int v1 = vect[1];
    alpha1 = e1p * Math::pgcd ( a1p, omegas );
    beta1 = -phi * v1;
    U1 = H * Vector3D ( -phi * u1, omegap * Math::pgcd ( a1p, omegas ), 0 );
  }
  
  {
    int i1p = i1 / Math::pgcd ( m_omega, i1 );
    int omegap = m_omega / Math::pgcd ( m_omega, i1 );
    int e1p = e1 / Math::pgcd ( Math::pgcd ( e1, m_omega ), f1 * omegap );
    int phi = f1 * omegap / Math::pgcd ( Math::pgcd ( e1, m_omega ), f1 * omegap );
    int omegas = m_omega / Math::pgcd ( Math::pgcd ( e1, m_omega ), f1 * omegap );
    Math::couple vect = Math::extendedEuclidean ( e1p, omegas );
    int u1 = vect[0];
    int v1 = vect[1];
    int psi = c1 * omegap * Math::pgcd ( e1p, omegas ) - b1 * phi * u1;
    int a1p = a1 / Math::pgcd ( Math::pgcd ( a1, psi ), Math::pgcd ( m_omega, omegas * b1 / Math::pgcd ( e1p, omegas ) ) );
    int psip = psi / Math::pgcd ( Math::pgcd ( a1, psi ), Math::pgcd ( m_omega, omegas * b1 / Math::pgcd ( e1p, omegas ) ) );
    int omegat = m_omega / Math::pgcd ( Math::pgcd ( a1, psi ), Math::pgcd ( m_omega, omegas * b1 / Math::pgcd ( e1p, omegas ) ) );
    int khi = omegas * b1 / Math::pgcd ( e1p, omegas ) / Math::pgcd ( Math::pgcd ( a1, psi ), Math::pgcd ( m_omega, omegas * b1 / Math::pgcd ( e1p, omegas ) ) );
    vect = Math::extendedEuclidean ( a1p, khi );
    int u2 = vect[0];
    int v2 = vect[1];
    vect = Math::extendedEuclidean ( Math::pgcd ( a1p, khi ), omegat );
    u2 *= vect[0];
    v2 *= vect[0];
    int w2 = vect[1];
    int k = -psip * v2;
    beta2 = -psip * w2;
    int alphas = Math::pgcd ( Math::pgcd ( a1p, khi ), omegat );
    int alphap = alphas * Math::pgcd ( e1p, omegas );
    alpha2 = alphap * i1p;
    gamma2 = -phi * v1 * alphas - k * e1p / Math::pgcd ( e1p, omegas );
    U2 = H * Vector3D ( -psip * u2, -phi * u1 * alphas + k * omegas / Math::pgcd ( e1p, omegas ), alphap * omegap );
  }
}

/**
 * computes the image of a point
 *
 * @param p the initial point
 *
 * @return the image of the point
 */
const Vector3D QAT3D::calculate ( const Vector3D p ) const
{
  return ( m_matrix * p + m_vector ) / m_omega;
}
