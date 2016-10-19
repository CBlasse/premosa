//
//  Parameter.cpp
//  proj_temp_cpp
//
//  Created by Corinna Blasse on 20.08.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include "Parameter.h"

#include <iostream>
extern "C" {
#include "utilities.h"
}

namespace ProjectionMethod {
  
  ParaSet::ParaSet (int argc, char * argv[], char * specifications []) {
    
    Process_Arguments(argc,argv,specifications,1);
    
//    ProjectionMethod::ParaSet parameter;
    
    if (Is_Arg_Matched("-r")) {         // grid size
      radius = Get_Int_Arg("-r");
      if (radius %2 != 0) {           // the grid size has to be even to enable a subsequent refinement
        radius++;
      }
    } else {
      radius = 20;
    }
    
    if (Is_Arg_Matched("-d")) {         // distance parameter
      distance = Get_Int_Arg("-d");
    } else {
      distance = 2;
    }
    
    if (Is_Arg_Matched("-l")) {         // number of cell layers
      layer = Get_Int_Arg("-l");
    } else {
      layer = 1;
    }
    
    if (Is_Arg_Matched("-t")) {         // threshold for brigth pixels
      threshold = Get_Int_Arg("-t");
    } else {
      threshold = 50;
    }
    
    if (Is_Arg_Matched("-hmd")) {      // export downsampled height map
      printHeightMap = true;
    } else {
      printHeightMap = false;
    }
    
    if (Is_Arg_Matched("-hmr")) {      // export height map
      printRealHeightMap = true;
    } else {
      printRealHeightMap = false;
    }
    
    if (Is_Arg_Matched("-mi")) {      // if the signal spreads several layers, then one can project the maximal intensities of the detected surface layers and its neighboring ones
      maxInterpolation = true;
    } else {
      maxInterpolation = false;
    }
    
    if (Is_Arg_Matched("-v")) {      // verbose parameter
      verbose = true;
    } else {
      verbose = false;
    }
  };
  
  ParaSet::~ParaSet () {
  }
 
  ExtendedParaSet::ExtendedParaSet (int argc, char * argv[], char * specifications []) {
    
    Process_Arguments(argc,argv,specifications,1);
    
    if (Is_Arg_Matched("-r")) {         // grid size
      radius = Get_Int_Arg("-r");
      if (radius %2 != 0) {           // the grid size has to be even to enable a subsequent refinement
        radius++;
      }
    } else {
      radius = 20;
    }
    
    if (Is_Arg_Matched("-d1")) {         // distance parameter
      distance1 = Get_Int_Arg("-d1");
    } else {
      distance1 = 0;
    }
  
    if (Is_Arg_Matched("-d2")) {         // distance parameter
      distance2 = Get_Int_Arg("-d2");
    } else {
      distance2 = 0;
    }
    
    if (Is_Arg_Matched("-t")) {         // threshold for brigth pixels
      threshold = Get_Int_Arg("-t");
    } else {
      threshold = 50;
    }
    
    if (Is_Arg_Matched("-hmd")) {      // export downsampled height map
      printHeightMap = true;
    } else {
      printHeightMap = false;
    }
    
    if (Is_Arg_Matched("-hmr")) {      // export height map
      printRealHeightMap = true;
    } else {
      printRealHeightMap = false;
    }
    
    if (Is_Arg_Matched("-v")) {      // verbose parameter
      verbose = true;
    } else {
      verbose = false;
    }
  };
  
  ExtendedParaSet::~ExtendedParaSet () {
  }

  
}