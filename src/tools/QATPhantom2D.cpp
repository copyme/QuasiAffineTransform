/**
 * @file   QATPhantom2D.cpp
 * @author David Coeurjolly <david.coeurjolly@liris.cnrs.fr>
 * @date   Thu Jun 19 10:59:44 2008
 *
 * @brief  2D Phantom
 *
 *
 */
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

#include "../common/mathematic.hpp"
#include "../tclap/CmdLine.h"
#include "../2D/image.hpp"

using namespace std;
using namespace TCLAP;

//Global Varialbe
int SIZE;

double zoneplate ( double  i, double  j)
{
  return ( 256* ( 127.5 + 127.5*cos ( ( 1440.0/M_PI ) / ( 1.0 + 512.0/ ( sqrt ( 8* ( ( i-SIZE/2 ) * ( i-SIZE/2 ) + ( j-SIZE/2 ) * ( j-SIZE/2 ) ) ) ) ) ) ) );
}

int main ( int argc, char** argv )
{
  try
    {
      //CmdLine Parser generator
      CmdLine cmd ( "Phantom image construction in dimension 2", ' ', "0.1" );
      ValueArg<int> imageSize ( "s","size","image size",true,512,"size" );
      cmd.add ( imageSize );

      // Parse the argv array.
      cmd.parse ( argc, argv );

      SIZE = imageSize.getValue();

      Image output ( SIZE,SIZE );
      double grey;
      Color c;
      for ( unsigned int i=0; i < SIZE; i++ )
        for ( unsigned int j=0; j< SIZE ; j++ )
          {
            grey = (int)rint(zoneplate ( i,j ));

            c.blue ( grey );
            c.red ( grey );
            c.green ( grey );
            output.setColor ( i,j,c );
          }

      output.write ( "zoneplate.png" );
    }
  catch ( ArgException &e )  // catch any exceptions
    {
      std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }


  return 0;
}
