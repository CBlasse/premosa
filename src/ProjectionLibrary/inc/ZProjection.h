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

#include "RasterizedImage.h"
#include "Parameter.h"

namespace ProjectionMethod {
  
  Array * MaxProjection_8bit (RasterizedImage & img, ParaSet & imgPara);
  Array * MaxProjection_16bit (RasterizedImage & img, ParaSet & imgPara);
  Array * Substack_UsingMask_8bit(RasterizedImage & img, ParaSet & imgPara, Array *median, bool onlyTopPixel);
  Array * Substack_UsingMask_16bit(RasterizedImage & img, ParaSet & imgPara, Array *median, bool onlyTopPixel);
  Array * MaxIndices_8bit(RasterizedImage & img, ParaSet & imgPara, Array *maxim);
  Array * MaxIndices_16bit(RasterizedImage & img, ParaSet & imgPara, Array *maxim);

  Array * GetSDImg (Array * image, int radius);
}

#endif
