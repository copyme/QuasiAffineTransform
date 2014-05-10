/**
 * @file   paving.cpp
 * @author Valentin Blot <valentin@Valentin-E21ter>
 * @date   Thu Jun 19 12:36:18 2008
 * 
 * @brief  Pavings (sets of points)
 * 
 * 
 */

#include "paving.hpp"

/** 
 * returns the number of points in the paving
 * 
 * 
 * @return size of the paving
 */
unsigned int Paving::size() const
{
  return m_points.size();
}

/** 
 * adds a point to the paving
 * 
 * @param point point to add
 */
void Paving::addPoint(const Vector2D point)
{
  m_points.push_back(point);
}

/** 
 * get a point of the paving - read-only version
 * 
 * @param index index of the point to get
 * 
 * @return the point
 */
const Vector2D Paving::operator[](unsigned int index) const
{
  return m_points[index] ;
}

/** 
 * get a point of the paving - read-write version
 * 
 * @param index index of the point to get
 * 
 * @return the point
 */
Vector2D &Paving::operator[](unsigned int index)
{
  return m_points[index] ;
}
