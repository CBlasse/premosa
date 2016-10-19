//
//  ZPlaneSelection.h
//  proj_temp_cpp
//
//  Created by Corinna Blasse on 26.07.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef ZPLANESELECTION_H
#define ZPLANESELECTION_H

extern "C" {
  #include "array.h"
}


namespace ImageClass{
  class Image;
}
namespace Parameter{
  struct ImgParameter;
}


namespace ZPlane {
  
  Array * PlaneSelection_8bit (ImageClass::Image & img, Parameter::ImgParameter & imgPara, Array *maxim, Array *index);
  Array * PlaneSelection_16bit (ImageClass::Image & img, Parameter::ImgParameter & imgPara, Array *maxim, Array *index);
  Array * PlaneSelection_DoubleRadius (ImageClass::Image & img, Parameter::ImgParameter & imgPara, Array *maxim, Array *index);
    
  Array * SelectPlanes_Variance_8bit(ImageClass::Image & img, Parameter::ImgParameter & imgPara, Array * substack);
  Array * SelectPlanes_Variance_16bit(ImageClass::Image & img, Parameter::ImgParameter & imgPara, Array * substack);
  Array * SelectPlanes_Entropy(ImageClass::Image & img, Parameter::ImgParameter & imgPara, Array * substack);
  Array * SelectPlanes_LOG(ImageClass::Image & img, Parameter::ImgParameter & imgPara, Array * substack);
  Array * SelectPlanes_LOG_basal(ImageClass::Image & img, Parameter::ImgParameter & imgPara, Array * substack);
}


#endif
