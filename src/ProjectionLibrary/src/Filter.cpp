//
//  Filter.cpp
//
//  Created by Corinna Blasse on 07/26/12.
//

#include "Filter.h"

#include <algorithm>
#include <iostream>
#include <vector>

extern "C" {
#include "histogram.h"
}

#include "HeightMapAnalysis.h"
#include "MylibValues.h"
#include "VectorFunctions.h"

#include "ShowArray.h"

using namespace MylibValues;
using namespace VectorFunctions;
using namespace ProjectionResults;


namespace ProjectionMethod {
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Standard median filter
   
   @param img:      input image
   @param imgPara:  struct with parameters
   @param original: current height map
   @param filterRadius:   Radius to define the 8-connected neighborhood
   */
  
  Array * MedianFilter(RasterizedImage & img, ParaSet & imgPara, Array *original, int filterRadius){
    
    // output image
    Array  *median = Make_Array_With_Shape(PLAIN_KIND,UINT32_TYPE,Coord2(original->dims[1], original->dims[0]));
    uint32 *med    = AUINT32(median);
    uint32 *orig    = AUINT32(original);
    
#ifdef DEVELOP
    Array  *variance = Make_Array_With_Shape(PLAIN_KIND,UINT32_TYPE,Coord2(original->dims[1], original->dims[0]));
    uint32 *var      = AUINT32(variance);
#endif
    
    Use_Reflective_Boundary();
    
    // frames defines the local neighborhood
    Frame * f = Make_Frame(original,Coord2(2*filterRadius+1,2*filterRadius+1),Coord2(filterRadius,filterRadius));
    Histogram * h = Make_Histogram(UVAL,img.GetDepth(),ValU(1),ValU(0));
    Place_Frame(f,0);
    
    for (Indx_Type p = 0; p < median->size; p++){
      Empty_Histogram(h);
      Histagain_Array(h,f,0);
      h->counts[orig[p]]--;              // excludes the height of p in the calculation of the median
      med[p] = Percentile2Bin(h,.5);    // median is given at 50th percentile
#ifdef DEVELOP
      var[p] = 10*Histogram_Sigma(h);
#endif
      Move_Frame_Forward(f);
    }
    
    Kill_Histogram(h);
    Kill_Frame(f);
    
#ifdef DEVELOP
    ShowArray(img, imgPara, "median.tif", median);
    ShowArray(img, imgPara, "variance.tif", variance);
    Free_Array(variance);
#endif
    
    if (imgPara.verbose){
      std::cout << "Planes filtered" << std::endl;
    }
    
    return (median);
  }
  
  
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Median filter that only modifies a pixel when it significantly differs from the median of its 8-connected neighborhood given by the filterRadius
   
   @param img:      input image
   @param imgPara:  struct with parameters
   @param original: current height map
   @param filterRadius:   radius to define the 8-connected neighborhood
   @param threshold: threshold to identify pixels that need to be modified
   */
  
  Array * LocalMedianFilter(RasterizedImage & img, ParaSet & imgPara, Array * original, int filterRadius, double threshold ){
    
    uint32 *orig    = AUINT32(original);
    Array  *median = Make_Array_With_Shape(PLAIN_KIND, UINT32_TYPE, Coord2(original->dims[1], original->dims[0]));
    uint32 *med    = AUINT32(median);
    
    Use_Extend_Boundary();
    
    Frame * f = Make_Frame(original,Coord2(2*filterRadius+1,2*filterRadius+1),Coord2(filterRadius, filterRadius));
    Histogram * h = Make_Histogram(UVAL,img.GetDepth(),ValU(1),ValU(0));
    
    
    Place_Frame(f,0);
    for (Indx_Type p = 0; p < median->size; p++){
      Empty_Histogram(h);
      Histagain_Array(h,f,0);
      h->counts[orig[p]]--;              // excludes the height of p in the calculation of the median
      
      if (abs(static_cast<int>(orig[p]) - Percentile2Bin(h,.5)) > threshold) {      // replaces a pixel only if the different between the pixel value and the median is greater than the threshold
        med[p] = Percentile2Bin(h,.5);
      } else {
        med[p] = orig[p];
      }
      Move_Frame_Forward(f);
    }
    
    Kill_Histogram(h);
    Kill_Frame(f);
    
    if (imgPara.verbose){
      std::cout << "Planes filtered" << std::endl;
    }
    
    return median;
  }
  
  
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Median filter that only modifies a pixel when it significantly differs from the median of its 8-connected neighborhood that is weighted by the variance
   
   @param img:        input image
   @param imgPara:    struct with parameters
   @param original:   current height map
   @param filterRadius:     radius to define the 8-connected neighborhood
   @param threshold:  threshold to identify pixels that need to be modified
   */
  Array * LocalMedianFilter_InclVar(RasterizedImage & img, ParaSet & imgPara, Array * original, int filterRadius, double threshold){
    
    uint32 *orig    = AUINT32(original);
    Array  *median = Make_Array_With_Shape(PLAIN_KIND, UINT32_TYPE, Coord2(original->dims[1], original->dims[0]));
    uint32 *med    = AUINT32(median);
    
    Use_Extend_Boundary();
    
    Frame * f = Make_Frame(original,Coord2(2*filterRadius+1,2*filterRadius+1),Coord2(filterRadius, filterRadius));
    
    // Frame to compute the variance in intensity of the current grid element
    Frame * varF = Make_Frame(img.GetImage(),Coord3(1, imgPara.radius,imgPara.radius),Coord3(0, 0, 0));
    Histogram * varH = Make_Histogram(UVAL,256,ValU(1),ValU(0));
    
    uint32 * data;
    
    Place_Frame(f,0);
    Place_Frame(varF, 0);
    
    
    for (Dimn_Type y = 0; y <= img.GetHeight()-imgPara.radius; y += imgPara.radius){
      for (Dimn_Type x = 0; x <= img.GetWidth()-imgPara.radius; x += imgPara.radius){
        
        Indx_Type pHM = (y/imgPara.radius) * original->dims[0] + (x/imgPara.radius);
        Place_Frame(f, pHM);
        data = (uint32 *) Frame_Values(f);
        std::vector<double> depths;
        double sum = 0;
        double n = 0;
        
        for (int i = 0; i < AForm_Size(f); i++) {
          Indx_Type p = Coord2IdxA(img.GetImage(), Coord3(data[i],y, x));
          Place_Frame(varF, p);
          Histagain_Array(varH,varF,0);
          
          sum += Histogram_Variance(varH) * data[i];
          n += Histogram_Variance(varH);
          for (int j=0; j< floor(Histogram_Variance(varH)/10); j++) {   // weight the current height with the variance
            depths.push_back(data[i]);
          }
        }
        double weightedMedian = VectorMedian(depths);
        
        if (abs(orig[pHM] - weightedMedian) > threshold ) {
          med[pHM] = floor(weightedMedian);
        } else {
          med[pHM] = orig[pHM];
        }
        
      }
    }
    
    Kill_Histogram(varH);
    Kill_Frame(varF);
    Kill_Frame(f);
    
    return median;
  }
  
  
  
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
  
  Array * LocalMedianFilter_Confidence(RasterizedImage & img, ParaSet & imgPara, Array * original, int filterRad, double threshold, Array * confidence){
    
    uint32 *orig    = AUINT32(original);
    Array  *median = Make_Array_With_Shape(PLAIN_KIND, UINT32_TYPE, Coord2(original->dims[1], original->dims[0]));
    uint32 *med    = AUINT32(median);
    
    Use_Extend_Boundary();
    
    Frame * f = Make_Frame(original,Coord2(2*filterRad+1,2*filterRad+1),Coord2(filterRad, filterRad));
    Frame * confF = Make_Frame(confidence,Coord2(2*filterRad+1,2*filterRad+1),Coord2(filterRad, filterRad));
    
    Place_Frame(f,0);
    Place_Frame(confF, 0);
    
    uint8 * confVals;
    uint32 * data;
    
    for (Indx_Type p = 0; p < median->size; p++){
      std::vector<double> heightVals;
      
      data = (uint32 *) Frame_Values(f);
      confVals = (uint8 *) Frame_Values(confF);
      
      for (int i = 0; i < AForm_Size(f); i++) {
        for (int j=0; j<confVals[i]; j++) {
          heightVals.push_back(data[i]);    // weight the height value with the 'likelihood' from the confidence image
        }
      }
      if (heightVals.size() == 0) {
        med[p] = orig[p];
      } else {
        med[p] = VectorMedian(heightVals);
      }
      
      Move_Frame_Forward(f);
      Move_Frame_Forward(confF);
    }
    
    Kill_Frame(f);
    Kill_Frame(confF);
    
    return median;
  }
  
  
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
  
  Array * LocalMedianFilter_3DConfidence(RasterizedImage & img, ParaSet & imgPara, Array * original, int filterRad, double threshold, Array * confidence){
    
    uint32 *orig    = AUINT32(original);
    Array  *median = Make_Array_With_Shape(PLAIN_KIND, UINT32_TYPE, Coord2(original->dims[1], original->dims[0]));
    uint32 *med    = AUINT32(median);
    float32 * confVals = AFLOAT32(confidence);
    Use_Extend_Boundary();
    
    Frame * f = Make_Frame(original,Coord2(2*filterRad+1,2*filterRad+1),Coord2(filterRad, filterRad));
//    Frame * confF = Make_Frame(confidence,Coord3(1, 2*filterRad+1,2*filterRad+1),Coord3(0, filterRad, filterRad));
    
    Place_Frame(f,0);
//    Place_Frame(confF, 0);
    
//    uint8 * confVals;
    uint32 * data;
    
    for (Indx_Type p = 0; p < median->size; p++){
      std::vector<double> heightVals;
      
      data = (uint32 *) Frame_Values(f);
//      confVals = (uint8 *) Frame_Values(confF);
      
      for (int i = 0; i < AForm_Size(f); i++) {
        
        Indx_Type q = data[i]* original->dims[1]*original->dims[0] + p;
        
        for (int j=0; j< confVals[q]; j++) {
          heightVals.push_back(data[i]);    // weight the height value with the 'likelihood' from the confidence image
        }
      }
      if (heightVals.size() == 0) {
        med[p] = orig[p];
      } else {
        med[p] = VectorMedian(heightVals);
      }
      
      Move_Frame_Forward(f);
//      Move_Frame_Forward(confF);
    }
    
    Kill_Frame(f);
//    Kill_Frame(confF);
    
    return median;
  }
  

  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Iterative filtering of the image using a median filter
   
   @param img:          input image
   @param imgPara:      struct with parameters
   @param original:     current height map
   @param filterRadii:  array with the filter radii to define the 8-connected neighborhood
   @param threshold:    threshold to identify pixels that need to be modified
   */
  
  /***********************
   * Iteration of the median filter using different radii
   **********************/
  
  Array * Filtering (RasterizedImage & img, ParaSet & imgPara, Array * original, int filterRadii [], double threshold) {
    Array * median = Copy_Array(original);
//    ComputeSmoothness(median);

    // iterate through the list of radii and filter the image
    for (int i=0; i < sizeof(filterRadii)/sizeof(*filterRadii); i++) {
      Free_Array(median);
      median = LocalMedianFilter(img, imgPara, original, filterRadii[i], threshold);
      
//      ComputeSmoothness(median);
      Free_Array(original);
      original = Copy_Array(median);
    }
    
#ifdef DEVELOP
    ShowArray(img, imgPara, "median.tif",median);
#endif
    if (imgPara.verbose){
      std::cout << "Planes filtered" << std::endl;
    }
    
    return median;
  }
  
  
  
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
  
  Array * Filtering_Confidence (RasterizedImage & img, ParaSet & imgPara, Array * original, int filterRadii [], double threshold, Array * confidence) {
    Array * median = Copy_Array(original);
    
    // iterate through the list of radii and filter the image
    for (int i=0; i < sizeof(filterRadii)/sizeof(*filterRadii); i++) {
      Free_Array(median);
      median = LocalMedianFilter_Confidence(img, imgPara, original, filterRadii[i], threshold, confidence);
      Free_Array(original);
      original = Copy_Array(median);
    }
    
#ifdef DEVELOP
    ShowArray(img, imgPara, "median.tif",median);
#endif
    if (imgPara.verbose){
      std::cout << "Planes filtered" << std::endl;
    }
    
    return median;
  }
  
  
  
}