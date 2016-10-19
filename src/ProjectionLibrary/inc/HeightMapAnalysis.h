//
//  HeightMapAnalysis.h
//
//  Created by Corinna Blasse on 08/30/12.
//

#ifndef HEIGHTMAPANALYSIS_H
#define HEIGHTMAPANALYSIS_H

extern "C" {
  #include "array.h"
}

namespace ProjectionResults {

  ////////////////////////////////////////////////////////////////////////////////
  /*
   Struct presenting the range of the height map
   */
  struct HeightMapRange {
    int minLayer;
    int maxLayer;
  };
  
  
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Computation of the range (minimal and maximal z-section) of the height map

   @param heightMap:  height map
   */
  HeightMapRange GetLevelRange(Array * heightMap);
  
  
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Computation of the overall smoothness of the height map
   (The smoothness is defined as the pixel-wise difference between the height map value and the median value of its neighbors)
   
   @param heightMap:  height map
   */
  double ComputeSmoothness(Array * heightMap);
  
}

#endif
