//
//  ZPlaneSelection.h
//  proj_temp_cpp
//
//  Created by Corinna Blasse on 26.07.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef ZPLANESELECTION_H
#define ZPLANESELECTION_H

#include <vector>

extern "C" {
  #include "array.h"
}

#include "RasterizedImage.h"
#include "Parameter.h"


namespace ProjectionMethod {
  
  struct Levels {
    std::vector<Array *> levels;
    Array * confidence;
  };
  
  std::vector<Array *> PlaneSelection_SeveralLevels_8bit (RasterizedImage & img, ParaSet & imgPara, Array *maxim, Array *index);
  std::vector<Array *> PlaneSelection_SeveralLevels_16bit (RasterizedImage & img, ParaSet & imgPara, Array *maxim, Array *index);
  Levels PlaneSelection_2Levels_16bit (RasterizedImage & img, ParaSet & imgPara, Array *maxim, Array *index);

  Array * PlaneSelection_8bit (RasterizedImage & img, ParaSet & imgPara, Array *maxim, Array *index);
  Array * PlaneSelection_16bit (RasterizedImage & img, ParaSet & imgPara, Array *maxim, Array *index);

  Array * PlaneSelection_DoubleRadius (RasterizedImage & img, ParaSet & imgPara, Array *maxim, Array *index);
  
  Array * SelectPlanes_Variance_8bit(RasterizedImage & img, ParaSet & imgPara, Array * substack);
  Array * SelectPlanes_Variance_16bit(RasterizedImage & img, ParaSet & imgPara, Array * substack);
  
  Array * SelectPlanes_Entropy(RasterizedImage & img, ParaSet & imgPara, Array * substack);
  Array * SelectPlanes_LOG(RasterizedImage & img, ParaSet & imgPara, Array * substack);
  Array * SelectPlanes_LOG_basal(RasterizedImage & img, ParaSet & imgPara, Array * substack);
  
  Array * BrightPixelDistribution (RasterizedImage & img, ParaSet & imgPara, Array *maxim, Array *index);
  
  
  Array * SelectPlanes_Variance_UsingHeightMap(RasterizedImage & img, ParaSet & imgPara, ExtendedParaSet & imgPara2, Array * substack);
  
}


#endif
