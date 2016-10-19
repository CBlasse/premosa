//
//  ZProjection.hsta
//  proj_temp_cpp
//
//  Created by Corinna Blasse on 25.07.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef ZPROJECTION_H
#define ZPROJECTION_H


extern "C" {
  #include "array.h"
}


namespace ImageClass{
  class Image;
}
namespace Parameter{
  struct ImgParameter;
}



namespace ZProjection {    
  
  Array * MaxProjection_8bit (ImageClass::Image & img, Parameter::ImgParameter & imgPara);
  Array * MaxProjection_16bit (ImageClass::Image & img, Parameter::ImgParameter & imgPara);

  Array * Substack_UsingMask_8bit(ImageClass::Image & img, Parameter::ImgParameter & imgPara, Array *median, bool onlyTopPixel);
  Array * Substack_UsingMask_16bit(ImageClass::Image & img, Parameter::ImgParameter & imgPara, Array *median, bool onlyTopPixel);
    
  Array * MaxIndices_8bit(ImageClass::Image & img, Parameter::ImgParameter & imgPara, Array *maxim);
  Array * MaxIndices_16bit(ImageClass::Image & img, Parameter::ImgParameter & imgPara, Array *maxim);
    }

#endif
