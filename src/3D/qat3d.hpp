#ifndef _QAT_H_
#define	_QAT_H_

#include <iostream>
#include <vector>
#include <math.h>

#include "matrix3x3.hpp"
#include "../common/mathematic.hpp"
#include "paving3d.hpp"
#include "image3d.hpp"

enum InterpolationType
{
  NO_BM,
  LINEAR ,
  NN
};


class QAT3D
{
 private:
  Matrix3x3 m_matrix;
  int m_omega;
  Vector3D m_vector;

  int a0, b0, c0, d0, e0, f0, g0, h0, i0, j0, k0, l0;
  int a1, b1, c1, d1, e1, f1, g1, h1, i1;
  int f10, f11, b10, b11;
  int p, p0, p1;
  int q, q0, q1;
  Matrix3x3 H, H0, H1, H2;
  int alpha0, alpha1, beta1, alpha2, beta2, gamma2;
  Vector3D U0, U1, U2;
  std::vector<Paving3D> pavings;

 public:
  QAT3D();
  QAT3D(Matrix3x3, int, Vector3D);

  bool isContracting() const;

  void inverse();

  const Paving3D determinePaving(int, int, int) const;
  const Paving3D determinePavingNoMultiply(int, int, int) const;

  const Paving3D determinePavingByBoundingRect(int, int, int) const;

  void setPaving(const Vector3D, const Vector3D);
  void setPaving(int, int, int, const Paving3D);
  const Paving3D getPaving(int, int, int) const;
  void determinePavings();
  void determinePavingsNoMultiply();
  void determinePavingsByBoundingRect();

  const Color backwardColorLinear(const Vector3D, Image3D) const;
  const Color backwardColorNN(const Vector3D, Image3D) const;

  
  const Image3D applyToImage(const Image3D, InterpolationType, bool, bool, bool, bool, bool);

  const Matrix3x3 getImageBound(const Image3D);

  void compute();

  const Vector3D calculate(const Vector3D) const;
};

#endif //_QAT_H_
