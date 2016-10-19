//
//  LocalIntegrals.h
//  CMakeProject
//
//  Created by blasse on 1/11/13.
//
//

#ifndef LOCALINTEGRALS_H
#define LOCALINTEGRALS_H

#include <iostream>

extern "C" {
#include "image.h"
#include "array.h"
#include "histogram.h"
#include "utilities.h"
#include "draw.h"
}

#include "Develop.h"
#include "HelperFunctions.h"



namespace Locally {
  Array * GetSDImg (Array * image);
  Array * GetAvgImg (Array * image);
  
  Array * LocalMedianFilter (Array * image);
}


#endif
