#ifndef _IMAGE_H_
#define _IMAGE_H_

#include "paving.hpp"
#include "../common/color.hpp"

#include <Magick++.h>
#include <string>

class Image
{
 private:
  Magick::Image img;
  bool wr;
  const Magick::PixelPacket *cacheRead;
  Magick::PixelPacket *cacheWrite;
  int w, h;

 public:
  Image();
  Image(int, int);
  Image(std::string);
  const Color getColor(int, int);
  const Color getFakeColor(int, int);
  void setColor(int, int, const Color);
  int height() const;
  int width() const;
  void write(const std::string);
  const Color colorOfPaving(const Paving&, const Vector2D&);
  void colorizePaving(const Paving&,  const Vector2D&, const Color&);
  void colorizePavingRemainder(const Paving&, const Paving&, const Vector2D&, 
		      const Color&, 
		      const Color&,
		      const Color &,
		      const Color &,
		      const Color &,
		      const int);

};

#endif // _IMAGE_H_
