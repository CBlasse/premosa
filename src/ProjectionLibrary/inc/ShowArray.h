//
//  ShowArray.h
//  proj_temp_cpp
//
//  Created by Corinna Blasse on 25.07.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef SHOWARRAY_H
#define SHOWARRAY_H

#include <string>

extern "C" {
  #include "array.h"
  #include "image.h"
  #include "draw.h"
}

#include "MylibValues.h"
#include "RasterizedImage.h"
#include "Parameter.h"


using namespace MylibValues;
using namespace ProjectionMethod;


namespace ProjectionResults {
  
  static Color_Bundle CYAN   = { MAX_PIX, ValU(0), ValU(150), ValU(150) };
  static Color_Bundle YELLOW = { MAX_PIX, ValU(150), ValU(150), ValU(0) };
  static Color_Bundle PURPLE = { MAX_PIX, ValU(150), ValU(0), ValU(150) };
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Method that presents the image arr with a grid overlap
   
   @param image:   input image stack
   @param imgPara: input parameter
   @param name:    name of the output file
   @param arr:     image which should be plotted
   */
  void ShowArray(RasterizedImage & image, ParaSet & imgPara, char * name, Array * arr);

}


#endif
