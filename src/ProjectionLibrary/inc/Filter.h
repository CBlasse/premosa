//
//  Filter.h
//
//  Created by Corinna Blasse on 07/20/15.
//
//

#ifndef FILTER_H
#define FILTER_H

extern "C" {
#include "array.h"
}

#include "RasterizedImage.h"
#include "Parameter.h"

namespace ProjectionMethod {
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Standard median filter
   
   @param img:      input image
   @param imgPara:  struct with parameters
   @param original: current height map
   @param filterRadius:   radius to define the 8-connected neighborhood
  */
  Array * MedianFilter(RasterizedImage & img, ParaSet & imgPara, Array *original, int filterRadius);
  
  
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Median filter that only modifies a pixel when it significantly differs from the median of its 8-connected neighborhood given by the filterRadius
    
   @param img:      input image
   @param imgPara:  struct with parameters
   @param original: current height map
   @param filterRadius:   radius to define the 8-connected neighborhood
   @param threshold: threshold to identify pixels that need to be modified
   */
  Array * LocalMedianFilter(RasterizedImage & img, ParaSet & imgPara, Array * original, int filterRadius, double threshold);
  
  
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Median filter that only modifies a pixel when it significantly differs from the median of its 8-connected neighborhood that is weighted by the variance
   
   @param img:        input image
   @param imgPara:    struct with parameters
   @param original:   current height map
   @param filterRadius:     radius to define the 8-connected neighborhood
   @param threshold:  threshold to identify pixels that need to be modified
  */
  Array * LocalMedianFilter_InclVar(RasterizedImage & img, ParaSet & imgPara, Array * original, int filterRadius, double threshold);

  
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Median filter that only modifies a pixel when it significantly differs from the median of its 8-connected neighborhood that is weighted by a likelihood
   
   @param img:        input image
   @param imgPara:    struct with parameters
   @param original:   current height map
   @param filterRadius:     radius to define the 8-connected neighborhood
   @param threshold:  threshold to identify pixels that need to be modified
   @param confidence: image indicating for each image a 'likelihood' that this represents the correct surface

   */
  Array * LocalMedianFilter_InclVar(RasterizedImage & img, ParaSet & imgPara, Array * original, int filterRadius, double threshold, Array * confidence);
  

  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Iterative filtering of the image using a median filter
   
   
   @param img:          input image
   @param imgPara:      struct with parameters
   @param original:     current height map
   @param filterRadii:  array with the filter radii to define the 8-connected neighborhood
   @param threshold:    threshold to identify pixels that need to be modified
   */
  Array * Filtering (RasterizedImage & img, ParaSet & imgPara, Array * original, int filterRadii [], double threshold);
  
  
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Iterative filtering of the image using a median filter, which integrates given likelihoods
   
   @param img:          input image
   @param imgPara:      struct with parameters
   @param original:     current height map
   @param filterRadii:  array with the filter radii to define the 8-connected neighborhood
   @param threshold:    threshold to identify pixels that need to be modified
   @param confidence:   image indicating for each image a 'likelihood' that this represents the correct surface
   */
  Array * Filtering_Confidence (RasterizedImage & img, ParaSet & imgPara, Array * original, int filterRadii [], double threshold, Array * confidence);
}



#endif
