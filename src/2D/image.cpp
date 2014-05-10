/**
 * @file   image.cpp
 * @author Valentin Blot <valentin@Valentin-E21ter>
 * @date   Thu Jun 19 13:44:38 2008
 *
 * @brief  images management (through libmagick++)
 *
 *
 */


#include "image.hpp"
#include <string>
#include <Magick++.h>

/**
 * default constructor
 *
 */
Image::Image()
{
}

/**
 * construct from size given by 2 integers
 *
 * @param width width of the image
 * @param height height of the image
 */
Image::Image(int width, int height)
{
    w = width;
    h = height;
    img = Magick::Image(Magick::Geometry(w, h), Magick::Color(0, 0, 0));
    img.type(Magick::TrueColorType);
    img.modifyImage();
    wr = true;
    cacheWrite = img.setPixels(0, 0, w, h);
}

/**
 * construct from file name
 *
 * @param fileName file name
 */
Image::Image(std::string fileName)
{
    img = Magick::Image(fileName);
    w = img.size().width();
    h = img.size().height();
    img.type(Magick::TrueColorType);
    img.modifyImage();
    wr = false;
    cacheRead = img.getConstPixels(0, 0, w, h);
}

/**
 * get color of a pixel
 *
 * @param i X-axis coordinate of the pixel
 * @param j Y-axis coordinate of the pixel
 *
 * @return color of the pixel
 */
const Color Image::getColor(int i, int j)
{
    if (wr)
    {
        img.syncPixels();
        cacheRead = img.getConstPixels(0, 0, w, h);
        wr = false;
    }
    const Magick::PixelPacket *pixel = cacheRead + w * (h - 1 - j) + i;
    return Color(pixel->red, pixel->green, pixel->blue);
}

/**
 * Create a fake color
 *
 * @param i X-axis coordinate of the pixel
 * @param j Y-axis coordinate of the pixel
 *
 * @return color of the pixel
 */
const Color Image::getFakeColor(int i, int j)
{
    int r,v,b;
    r = (i%2)* 65535;
    v = (j%2)* 65535;
    b = ((i+j)%2)* 65535;
    return Color(r,v,b);
}


/**
 * set color of a pixel
 *
 * @param i X-axis coordinate of the pixel
 * @param j Y-axis coordinate of the pixel
 * @param color color to give to the pixel
 */
void Image::setColor(int i, int j, const Color color)
{
    if (!wr)
    {
        cacheWrite = img.setPixels(0, 0, w, h);
        wr = true;
    }
    Magick::PixelPacket *pixel = cacheWrite + w * (h - 1 - j) + i;
    pixel->red = color.red();
    pixel->green = color.green();
    pixel->blue = color.blue();
}

/**
 * returns the height of the image
 *
 *
 * @return height of the image
 */
int Image::height() const
{
    return h;
}

/**
 * returns the width of the image
 *
 *
 * @return width of the image
 */
int Image::width() const
{
    return w;
}

/**
 * writes the image to a file
 *
 * @param fileName file name of the file to put the image in
 */
void Image::write(const std::string fileName)
{
    if (wr)
    {
        img.syncPixels();
        img.write(fileName);
    }
}

/**
 * gets the average color of a translated paving in the image
 *
 * @param P the paving
 * @param vect the translating vector
 *
 * @return the average color
 */
const Color Image::colorOfPaving(const Paving &P, const Vector2D &vect)
{
    int count = 0;
    Color color(0, 0, 0);
    for (unsigned k = 0; k < P.size(); k++)
    {
        Vector2D P2 = P[k] + vect;
        if (0 <= P2[0] && P2[0] < w && 0 <= P2[1] && P2[1] < h)
        {
            color += getColor(P2[0], P2[1]);
            count++;
        }
    }
    if (count)
        color /= count;
    return color;
}

#include<iostream>

/**
 * sets the color of pixels belonging to a translated paving
 *
 * @param P the paving
 * @param vect the translating vector
 * @param color the color to give to the paving 
 */
void Image::colorizePaving(const Paving &P, const Vector2D &vect, 
			   const Color &color)
{
  
  for (unsigned k = 0; k < P.size(); k++)
    {
	Vector2D P2 = P[k] + vect;
        if (0 <= P2[0] && P2[0] < w && 0 <= P2[1] && P2[1] < h)
	  setColor(P2[0], P2[1], color);
    }
}


/**
 * sets the color of pixels belonging to a translated paving
 *
 * @param P the paving
 * @param vect the translating vector
 * @param color the color to give to the paving 
 */
void Image::colorizePavingRemainder(const Paving &P, const Paving &rP, const Vector2D &vect, 
			   const Color &color,
			   const Color &cN,
			   const Color &cS,
			   const Color &cE,
			   const Color &cO,
			   const int omega)
{
    ///FIXME on trouve les max et on normalise par rapport au reste X
  
  int maxX=0, maxY=0;
    
  std::cout <<"rP = { ";
  for ( unsigned int k = 0 ; k< rP.size(); k++)
    {
      std::cout <<"("<<rP[k][0] << ","<<rP[k][1]<<") ";
      if (rP[k][0]> maxX)
	maxX = rP[k][0];
      if (rP[k][1]> maxY)
	maxY = rP[k][1];

    }
  std::cout << " }"<<std::endl;

  double alpha;
  
  for (unsigned k = 0; k < P.size(); k++)
    {
	Vector2D P2 = P[k] + vect;

	alpha = (double)(rP[k][1])  / omega;
	
        if (0 <= P2[0] && P2[0] < w && 0 <= P2[1] && P2[1] < h)
            setColor(P2[0], P2[1], color*alpha);
    }
}
