Introduction
============

This package implements fast image transformations in dimension 2 and 3 with quasi-affine functions.

The code follows the following papers (https://liris.cnrs.fr/david.coeurjolly/publications.html):

- Quasi-Affine Transformation in 3-D: Theory and Algorithms , D. Coeurjolly, V. Blot, M.A. Jacob-DaCol, 13th International Workshop on Combinatorial Image Analysis, Mexico, Springer Verlag, LNCS, 2009.
- Quasi-Affine Transformation in Higher Dimension , V. Blot and D. Coeurjolly, 15th IAPR International Workshop on Discrete Geometry for Computer Imagery, Montreal, Springer LNCS,2009.





LICENSE
=======

- GPLv2, see LICENSE file


TODO
====

Implement  usePeriodicity and noMultiply in 3D

Files
=====

Arithmetic Kernel:
common/mathematic.cpp
common/mathematic.hpp

Main classes (data-structure) :
common/color.cpp
common/color.hpp

2D
--

Main classes:
2D/image.cpp
2D/image.hpp
2D/paving.cpp
2D/paving.hpp
2D/matrix2x2.cpp
2D/matrix2x2.hpp
2D/vector2d.cpp
2D/vector2d.hpp

Transformation algorithm:
2D/qat.cpp
2D/qat.hpp

Main : 
2D/qat-2D.cpp

3D
--

Main classes:
3D/image3d.cpp
3D/image3d.hpp
3D/paving3d.cpp
3D/paving3d.hpp
3D/matrix3x3.cpp
3D/matrix3x3.hpp
3D/vector3d.cpp
3D/vector3d.hpp

Transformation algorithm:
3D/qat3d.cpp
3D/qat3d.hpp

Main : 
3D/qat-3D.cpp
