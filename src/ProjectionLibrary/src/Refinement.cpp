//
//  Refinement.cpp
//  proj_temp_cpp
//
//  Created by Corinna Blasse on 16.08.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>

#include "Refinement.h"


namespace ProjectionMethod {
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Decomposition of the height map using the quad-tree principle (each  pixel is split into four)
   
   @param heightMap:  height map
   */
  
  Array * DecomposeImage (Array * original){
    Dimn_Type x, y;
    Indx_Type q;
    
    uint32 * val = AUINT32(original);
    uint32 value;
    
    Array * refImg = Make_Array_With_Shape(PLAIN_KIND, UINT32_TYPE, Coord2(2*original->dims[1], 2*original->dims[0]));
    uint32 * refVal = AUINT32(refImg);
    
    for (y=0; y < original->dims[1]; y++) {
      for (x=0; x < original->dims[0]; x++) {
        value = val[Coord2IdxA(original, Coord2(y, x))];
        
        q = Coord2IdxA(refImg, Coord2(2*y, 2*x));
        refVal[q] = value;
        
        q = Coord2IdxA(refImg, Coord2(2*y+1, 2*x));
        refVal[q] = value;
        
        q = Coord2IdxA(refImg, Coord2(2*y, 2*x+1));
        refVal[q] = value;
        
        q = Coord2IdxA(refImg, Coord2(2*y+1, 2*x+1));
        refVal[q] = value;
      }
    }
    
    return refImg;

  }
}