#ifndef _IMAGE3D_H_
#define _IMAGE3D_H_

#include "paving3d.hpp"
#include "../common/color.hpp"

#include <string>
#include <vol.h>

class Image3D
{
 private:
  Vol vol;
	
 public:
  Image3D();
  Image3D(int, int, int);
  Image3D(std::string);
  const Color getColor(int, int, int) const;
  const Color getFakeColor(int, int, int) const;
  void setColor(int, int, int, const Color);
  int height() const;
  int width() const;
  int depth() const;
  void write(const std::string);
  const Color colorOfPaving(const Paving3D, const Vector3D) const;
  void colorizePaving(const Paving3D, const Vector3D, const Color);
	
};


/**
 * get color of a pixel
 *
 * @param i X-axis coordinate of the pixel
 * @param j Y-axis coordinate of the pixel
 * @param k Z-axis coordinate of the pixel
 *
 * @return color of the pixel
 */
inline const Color Image3D::getColor ( int i, int j, int k ) const
  {
    if ( i < 0 || i >= width() || j < 0 || j >= height() || k < 0 || k >= depth() )
      return Color ( 0, 0, 0 );
    voxel col = vol ( vol.minX() + i, vol.minY() + j, vol.minZ() + k );
    return Color ( col, col, col );
  }

/**
 * get fake color of a pixel (except for background voxels)
 *
 * @param i X-axis coordinate of the pixel
 * @param j Y-axis coordinate of the pixel
 * @param k Z-axis coordinate of the pixel
 *
 * @return fake color of the pixel
 */
inline const Color Image3D::getFakeColor ( int i , int j , int k ) const
  {
    int r,v,b;
    r = ( i%2 ) * 65535;
    v = ( j%2 ) * 65535;
    b = ( ( i+j+k ) %2 ) * 65535;
    voxel col = vol ( vol.minX() + i, vol.minY() + j, vol.minZ() + k );
    if ( col == 0 )
      return Color ( col,col,col );
    else
      return Color ( r,v,b );
  }


/**
 * set color of a pixel
 *
 * @param i X-axis coordinate of the pixel
 * @param j Y-axis coordinate of the pixel
 * @param k Z-axis coordinate of the pixel
 * @param color color to give to the pixel
 */
inline void Image3D::setColor ( int i, int j, int k, const Color color )
{
  voxel col = ( color.red() + color.green() + color.blue() ) / 3;
  vol ( vol.minX() + i, vol.minY() + j, vol.minZ() + k ) = col;
}

/**
 * returns the width of the image
 *
 *
 * @return width of the image
 */
inline int Image3D::width() const
  {
    return vol.sizeX();
  }

/**
 * returns the height of the image
 *
 *
 * @return height of the image
 */
inline int Image3D::height() const
  {
    return vol.sizeY();
  }

/**
 * returns the depth of the image
 *
 *
 * @return depth of the image
 */
inline int Image3D::depth() const
  {
    return vol.sizeZ();
  }



#endif // _IMAGE3D_H_
