/**
 * @file   image.cpp
 * @author Valentin Blot <valentin@Valentin-E21ter>
 * @date   Thu Jun 19 13:44:38 2008
 *
 * @brief  images management (through libmagick++)
 *
 *
 */


#include "image3d.hpp"

/**
 * default constructor
 *
 */
Image3D::Image3D()
{
}

/**
 * construct from size given by 2 integers
 *
 * @param width width of the image
 * @param height height of the image
 * @param depth depth of the image
 */
Image3D::Image3D ( int width, int height, int depth )
{
  vol = Vol ( width, height, depth, 0 );
}

/**
 * construct from file name
 *
 * @param fileName file name
 */
Image3D::Image3D ( std::string fileName )
{
  vol = Vol ( fileName.c_str() );
}

/**
 * writes the image to a file
 *
 * @param fileName file name of the file to put the image in
 */
 void Image3D::write ( const std::string fileName )
{
  vol.dumpVol ( fileName.c_str() );
}

/**
 * gets the average color of a translated paving in the image
 *
 * @param P the paving
 * @param vect the translating vector
 *
 * @return the average color
 */
const Color Image3D::colorOfPaving ( const Paving3D P, const Vector3D vect ) const
  {
    int count = 0;
    Color color ( 0, 0, 0 );
    for ( unsigned k = 0; k < P.size(); k++ )
      {
        Vector3D P2 = P[k] + vect;
        if ( 0 <= P2[0] && P2[0] < width() && 0 <= P2[1] && P2[1] < height() && 0 <= P2[2] && P2[2] < depth() )
          {
            color += getColor ( P2[0], P2[1], P2[2] );
            count++;
          }
      }
    if ( count )
      color /= count;
    return color;
  }

/**
 * sets the color of pixels belonging to a translated paving
 *
 * @param P the paving
 * @param vect the translating vector
 * @param color the color to give to the paving
 */
void Image3D::colorizePaving ( const Paving3D P, const Vector3D vect, const Color color )
{
  for ( unsigned k = 0; k < P.size(); k++ )
    {
      Vector3D P2 = P[k] + vect;
      if ( 0 <= P2[0] && P2[0] < width() && 0 <= P2[1] && P2[1] < height() && 0 <= P2[2] && P2[2] < depth() )
        setColor ( P2[0], P2[1], P2[2], color );
    }
}
