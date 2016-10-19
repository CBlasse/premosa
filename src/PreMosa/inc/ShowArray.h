//
//  ShowArray.h
//  proj_temp_cpp
//
//  Created by Corinna Blasse on 25.07.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef SHOWARRAY_H
#define SHOWARRAY_H


extern "C" {
#include "array.h"
#include "draw.h"
}

#include "HelperFunctions.h"
using namespace HelperFunctions;

namespace ImageClass{
  class Image;
}
namespace Parameter{
  struct ImgParameter;
}


static Color_Bundle CYAN   = { MAX_PIX, ValU(0), ValU(150), ValU(150) };
static Color_Bundle YELLOW = { MAX_PIX, ValU(150), ValU(150), ValU(0) };
static Color_Bundle PURPLE = { MAX_PIX, ValU(150), ValU(0), ValU(150) };

void ShowArray(ImageClass::Image & img, Parameter::ImgParameter & imgPara, char * name, Array * arr);


#endif
