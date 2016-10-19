//
//  AdjustHistogram.h
//  CMakeProject
//
//  Created by blasse on 1/18/13.
//
//

#ifndef ADJUSTHISTOGRAM_H
#define ADJUSTHISTOGRAM_H

#include <iostream>

extern "C" {
#include "array.h"
#include "histogram.h"
}

#include "Parameter.h"


namespace Parameter {
  struct ImgParameter;
}


namespace AdjustHistogram {
  
  int Percentile2SmallestBin (Histogram * h, double percentile);
  double Bin2PercentileReversed (Histogram * h, int b);
  
  Array * HistogramFitting (Array * projection, Parameter::ImgParameter & imgPara, Array * reference);
  Array * SmoothLowIntensityRegions (Array * projection, Parameter::ImgParameter & imgPara);
  Array * ImproveBorderRegions (Array * image, Parameter::ImgParameter & imgPara, Array * reference);
  
}


#endif 
