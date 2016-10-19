//
//  Correlation.h
//  CMakeProject
//
//  Created by blasse on 4/4/13.
//
//

#ifndef CORRELATION_H
#define CORRELATION_H

#include <iostream>

extern "C" {
#include "array.h"
}


namespace Correlation {
  
  double CorrelationOf2Images (Array * img1, Array img2);
  
}

#endif 