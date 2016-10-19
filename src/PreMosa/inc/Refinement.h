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



namespace Refinement {
  Array * DecomposeImage (Array * original);
}


#endif
