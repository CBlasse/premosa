//
//  Interpolation.h
//  proj_temp_cpp
//
//  Created by Corinna Blasse on 27.07.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef INTERPOLATION_H
#define INTERPOLATION_H

extern "C" {
#include "array.h"
}

namespace ImageClass {
  class Image;
}
namespace Parameter {
  struct ImgParameter;
}


namespace Interpolation {
  Array * InterpolateCorners(ImageClass::Image & img, Parameter::ImgParameter & imgPara, Array *median);

  Array * InterpolatePlanes_8bit(ImageClass::Image & img,  Parameter::ImgParameter & imgPara, Array *maxim, Array *corners);
  Array * InterpolatePlanes_16bit(ImageClass::Image & img,  Parameter::ImgParameter & imgPara, Array *maxim, Array *corners);
}

#endif
