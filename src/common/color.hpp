#ifndef _COLOR_H_
#define _COLOR_H_

class Color
{
 private:
  int r;
  int g;
  int b;
	
 public:
  Color();
  Color(int, int, int);
  int red() const;
  void red(int);
  int green() const;
  void green(int);
  int blue() const;
  void blue(int);
  const Color operator*(const float) const;
  Color &operator+=(const Color);
  Color &operator/=(int);
  Color &operator/=(double);
	
};


/** 
 * default constructor
 * 
 */
inline Color::Color()
{
}

/** 
 * construct from 3 integers
 * 
 * @param red red value
 * @param green green value
 * @param blue blue value
 */
inline Color::Color(int red, int green, int blue)
{
  r = red;
  g = green;
  b = blue;
}

/** 
 * get red value
 * 
 * 
 * @return red value
 */
inline int Color::red() const
{
  return r;
}

/** 
 * set red value
 * 
 * @param val red value
 */
inline void Color::red(int val)
{
  r = val;
}

/** 
 * get green value
 * 
 * 
 * @return green value
 */
inline int Color::green() const
{
  return g;
}

/** 
 * set green value
 * 
 * @param val green value
 */
inline void Color::green(int val)
{
  g = val;
}

/** 
 * get blue value
 * 
 * 
 * @return blue value
 */
inline int Color::blue() const
{
  return b;
}

/** 
 * set blue value
 * 
 * @param val blue value
 */
inline void Color::blue(int val)
{
  b = val;
}

/** 
 * computes the product of current color by a floating point value
 * 
 * @param val the value
 * 
 * @return the result
 */
inline const Color Color::operator*(const float val) const
{
  return Color(r * val, g * val, b * val);
}

/** 
 * add another color to current color
 * 
 * @param color the other color
 * 
 * @return the modified color
 */
inline Color &Color::operator+=(const Color color)
{
  r += color.r;
  g += color.g;
  b += color.b;
  return *this;
}

/** 
 * divide the color by a value
 * 
 * @param val the value
 * 
 * @return the modified color
 */
inline Color &Color::operator/=(int val)
{
  r /= val;
  g /= val;
  b /= val;
  return *this;
}

/** 
 * divide the color by a floating point value
 * 
 * @param val the value
 * 
 * @return the modified color
 */
inline Color &Color::operator/=(double val)
{
  r /= val;
  g /= val;
  b /= val;
  return *this;
}


#endif // _COLOR_H_
