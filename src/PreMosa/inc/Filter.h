//
//  Filter.h
//  proj_temp_cpp
//
//  Created by Corinna Blasse on 26.07.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef FILTER_H
#define FILTER_H

extern "C" {
#include "array.h"
}


namespace ImageClass {
  class Image;
}
namespace Parameter {
  struct ImgParameter;
}


namespace Filter {
  
  Array * MedianFilter(ImageClass::Image & img, Parameter::ImgParameter & imgPara, Array *original, int depth, int rad);
  Array * LocalMedianFilter(ImageClass::Image & img, Parameter::ImgParameter & imgPara, Array * original, int depth, int filterRad);
  Array * LocalMedianFilter_InclVar(ImageClass::Image & img, Parameter::ImgParameter & imgPara, Array * original, int depth, int filterRad);
  Array * Filtering (ImageClass::Image & img, Parameter::ImgParameter & imgPara, Array * original, int depth, int filterSizes []);
  
}



#endif
