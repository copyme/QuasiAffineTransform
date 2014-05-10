/**
 * @file   qat-3D.cpp
 * @author Valentin Blot <valentin@Valentin-E21ter>
 * @author David Coeurjolly <david.coeurjolly@liris.cnrs.fr>
 * @date   Thu Jun 19 10:59:44 2008
 *
 * @brief  main file
 *
 *
 */
#include <iostream>
#include <fstream>

#include "../common/mathematic.hpp"
#include "../tclap/CmdLine.h"
#include "qat3d.hpp"

using namespace std;
using namespace TCLAP;



/**
 * Compute the PSNR between two images (used in the f.f^-1 test)
 *
 * @param im1 image 1
 * @param im2 image 2
 *
 * @return the PSNR value in dB
 */
double psnr ( Image3D &im1, Image3D &im2 )
{
  double a,b,diff = 0;
  Color c;
  double max=0;
  double nb=0;

  int Xi = (im2.width() - im1.width()) / 2;
  int Yi = (im2.height() - im1.height()) / 2;
  int Zi = (im2.depth() - im1.depth()) / 2;

  for ( unsigned int i=0 ; i < im1.width() ; i++ )
    for ( unsigned int j=0 ; j < im1.height() ; j++ )
      for ( unsigned int k=0 ; k < im1.depth() ; k++ )
        {
          c = im1.getColor ( i,j,k );
          a = c.green() + c.red() + c.blue();
          c = im2.getColor ( Xi + i, Yi + j, Zi + k );
          b = c.green() + c.red() + c.blue();

          if ( a>max ) max = a;
          diff += ( a -b ) * ( a-b );
          nb++;
        }
  return 20*log10 ( max/sqrt ( diff/nb ) );
}

/**
 * main procedure
 *
 * @param argc # of args
 * @param argv args
 *
 * @return no error
 */
int main ( int argc, char** argv )
{

  try
    {
      //CmdLine Parser generator
      CmdLine cmd ( "Quasi Affine Transform in dimension 3", ' ', "0.3" );
      vector<Arg *> xorlist;

      //Option (not required)
      SwitchArg backwardSwitch ( "l","linearbackward","Bilinear Backward Mapping", false );
      xorlist.push_back ( &backwardSwitch );
      SwitchArg backwardNNSwitch ( "n","NNbackward","Nearest Neighbor Backward Mapping", false );
      xorlist.push_back ( &backwardNNSwitch );
      SwitchArg naiveSwitch ( "","naive","Naive BoundingRect method",false );
      xorlist.push_back ( &naiveSwitch );
      SwitchArg periodicitySwitch ( "p","periodicity","Use paving periodicity", true );
      xorlist.push_back ( &periodicitySwitch );
      cmd.xorAdd ( xorlist );

      SwitchArg fakeColorSwitch ( "f","fake_color","Output fake colors to illustrate pavings (non contracting AQA only)", false );
      cmd.add ( fakeColorSwitch );
      SwitchArg nomultiplySwitch ( "m","no_multiplication","No multiplications in the paving computation", false );
      cmd.add ( nomultiplySwitch );
      SwitchArg compositionSwitch ( "","composition","Composition test: f.f^{-1}", false );
      cmd.add ( compositionSwitch );
      ValueArg<string> transformFile ( "t","transform","The transform file name (this file should contain in this order : omega, a, b, c, d, e, f, g, h, i, j, k, l, separated with spaces, where the quasi-affine transform is : (ax + by + cz + j, dx + ey + fz + k, gx + hy + iz + l) / omega)",true,"","file name" );
      cmd.add ( transformFile );
      ValueArg<string> outFile ( "o","output","The output image file name",true,"","file name" );
      cmd.add ( outFile );
      ValueArg<string> inFile ( "i","input","The input image file name",true,"","file name" );
      cmd.add ( inFile );
      // Parse the argv array.
      cmd.parse ( argc, argv );

      // Get the value parsed by each arg.
      // Get the value parsed by each arg.
      InterpolationType interp = NO_BM;
      if ( backwardSwitch.getValue() )
        interp = LINEAR;
      else
        if ( backwardNNSwitch.getValue() )
          interp = NN;

      bool useBoundingRect = naiveSwitch.getValue();
      bool noMultiply = nomultiplySwitch.getValue();
      bool usePeriodicity = periodicitySwitch.getValue();
      bool fakeColor = fakeColorSwitch.getValue();

      cout <<"Cmd Line Test: [ ";

      cout << "Backward mapping : "<< interp<<", ";

      cout << "BoundingRect : ";
      if ( useBoundingRect )
        cout << "Y, ";
      else
        cout << "N, ";

      cout << "noMultiply : ";
      if ( noMultiply )
        cout << "Y, ";
      else
        cout << "N, ";

      cout << "Periodicity : ";
      if ( usePeriodicity )
        cout << "Y";
      else
        cout << "N";
      cout << " ]"<<endl;


      // Aquisition de données
      int o, a, b, c, d, e, f, g, h, i, j, k, l;
      fstream transform ( transformFile.getValue().c_str(), fstream::in );
      transform >> o>> a >> b >> c >> d >> e >> f >> g >> h >> i >> j >> k >> l;
      transform.close();

      Image3D initialImage ( inFile.getValue() );

      QAT3D qat ( Matrix3x3 ( a, b, c, d, e, f, g, h, i ), o, Vector3D ( j, k, l ) );
      QAT3D qat2 ( Matrix3x3 ( a, b, c, d, e, f, g, h, i ), o, Vector3D ( j, k, l ) );

      if ( compositionSwitch.getValue() )
        {
          Image3D finalImage = qat.applyToImage ( initialImage, interp, useBoundingRect, usePeriodicity, noMultiply, fakeColor, false );

          cout << "Inverse computation..."<<endl;

          Image3D inverse = qat2.applyToImage ( finalImage, interp, useBoundingRect, usePeriodicity, noMultiply, fakeColor, true );
          inverse.write ( outFile.getValue() );

          double db = psnr ( initialImage, inverse );
          cout << "PSNR = "<<db<<" db"<<endl;
        }
      else
        {
          Image3D finalImage = qat.applyToImage ( initialImage, interp, useBoundingRect, usePeriodicity, noMultiply, fakeColor, false );
          finalImage.write ( outFile.getValue() );
        }


    }
  catch ( ArgException &e )  // catch any exceptions
    {
      std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }

  return 0;
}
