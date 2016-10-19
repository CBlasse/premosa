//
//  HeightMapAnalysis.cpp
//
//  Created by Corinna Blasse on 30/08/12.
//

#include "HeightMapAnalysis.h"


#include <cmath>
#include <iostream>


extern "C" {
#include "draw.h"
#include "histogram.h"
#include "image.h"
}

#include "MylibValues.h"

#include "Develop.h"

using namespace MylibValues;


namespace ProjectionResults {
  
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Computation of the range (minimal and maximal z-section) of the height map
   
   @param heightMap:  height map
   */
		
  HeightMapRange GetLevelRange(Array * heightMap){
    
    uint32 * hmVals = AUINT32(heightMap);
    HeightMapRange range;   // set inital values, which will be definitely replaced
    range.maxLayer = 0;
    range.minLayer = 1000000;
    
    for (int i=0; i<heightMap->size; i++) {
      if (hmVals[i] > range.maxLayer) {
        range.maxLayer = hmVals[i];
      } else if (hmVals[i] < range.minLayer){
        range.minLayer = hmVals[i];
      }
    }
    return range;
  }
  
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Computation of the overall smoothness of the height map
   (The smoothness is defined as the pixel-wise difference between the height map value and the median value of its neighbors)
   
   @param heightMap:  height map
   */
  
  double ComputeSmoothness(Array * heightMap) {
    
    Array  * distances = Make_Array_With_Shape(PLAIN_KIND,UINT32_TYPE,Coord2(heightMap->dims[1], heightMap->dims[0]));
    uint32 * distVals   = AUINT32(distances);
    uint32 * hmVals = AUINT32(heightMap);
    
    Use_Reflective_Boundary();
    
    HeightMapRange range = GetLevelRange(heightMap);
    
    // frame defines the local neighborhood
    Frame * f = Make_Frame(heightMap,Coord2(3,3),Coord2(1,1));
    Histogram * h = Make_Histogram(UVAL,range.maxLayer,ValU(1),ValU(0));
    Place_Frame(f,0);
    
    
    for (Indx_Type i=0; i< heightMap->size; i++) {
      
      Empty_Histogram(h);
      Histagain_Array(h, f, 0);
      
      int middleBin = Value2Bin(h, ValU(hmVals[i]));    // To determine the median value of the pixel's neighborhood, we exlude the current pixel i from the histogram
      h->counts[middleBin]--;
      
      distVals[i] = std::abs(static_cast<double>(hmVals[i]) - static_cast<double>(Percentile2Bin(h, 0.5)));
    }
    
#ifdef DEVELOP
    Write_Image("HeightMapSmoothness", distances, DONT_PRESS);
#endif
    
    Empty_Histogram(h);
    Histagain_Array(h, distances, 0);
    
    double meanDistance = Histogram_Mean(h);
    std::cout << "Smoothness:\n   Mean distance: " << meanDistance << "\n   Sd: " << Histogram_Sigma(h) << std::endl;
    
    Free_Histogram(h);
    Free_Frame(f);
    Free_Array(distances);
    
    return meanDistance;
  }
}

