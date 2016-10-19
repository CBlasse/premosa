//
//  Refinement.h
//  proj_temp_cpp
//
//  Created by Corinna Blasse on 16.08.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef REFINEMENT_H
#define REFINEMENT_H

extern "C" {
  #include "array.h"
}

namespace ProjectionMethod {
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Decomposition of the height map using the quad-tree principle (each  pixel is split into four)
   
   @param heightMap:  height map
   */
  Array * DecomposeImage (Array * heightMap);
}


#endif
