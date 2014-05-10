#ifndef _QAT_H_
#define	_QAT_H_

#include <iostream>
#include <vector>
#include <math.h>

#include "matrix2x2.hpp"
#include "../common/mathematic.hpp"
#include "paving.hpp"
#include "image.hpp"


enum InterpolationType
{
  NO_BM,
  LINEAR ,
  NN
};

class QAT
{
 private:
  Matrix2x2 m_matrix;
  int m_omega;
  Vector2D m_vector;

  int a0, b0, c0, d0, e0, f0;
  int a1, b1, c1, d1;
  int b10, b11;
  Matrix2x2 H, H0;
  int alpha0, alpha1, beta1;
  Vector2D U0, U1;
  std::vector<Paving> pavings;
  std::vector<Paving> pavingsRemainder; // remainders of the canonical paving points
  
 public:
  QAT();
  QAT(Matrix2x2, int, Vector2D);

  bool isContracting() const;
  bool isInversible() const;

  void inverse();

  const Paving determinePaving(int, int) const;
  const Paving determinePavingNoMultiply(int, int) const;
  const Paving determinePavingByBoundingRect(int, int) const;

  void setPaving(const Vector2D, const Vector2D);
  void setPavingRemainder ( const Vector2D I, const Vector2D Rem, const Vector2D P );
  
  void setPaving(int, int, const Paving);
  const Paving getPaving(int, int) const;
  const Paving getPavingRemainder ( int , int ) const;
  
  void determinePavings();
  void determinePavingsWithRemainders();
  void determinePavingsNoMultiply();
  void determinePavingsByBoundingRect();

  const Color backwardColorLinear(const Vector2D, Image) const;
  const Color backwardColorNN(const Vector2D, Image) const;
  
  const Image applyToImage(Image, InterpolationType , bool, bool, bool,bool, bool);

  const Matrix2x2 getImageBound(const Image);

  void compute();

  const Vector2D calculate(const Vector2D) const;
  const Vector2D calculateRemainder ( const Vector2D) const;
   
};

#endif //_QAT_H_
