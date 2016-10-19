//
//  MylibValues.h
//
//  Created by Corinna Blasse on 04/01/15.
//
//

#ifndef MYLIBVALUES_H
#define MYLIBVALUES_H

extern  "C" {
#include "mylib.h"
}

namespace MylibValues {
  
  Value ValU(uint64 value);
  Value ValU8(uint8 value);
  Value ValI(int64 value);
  Value ValF(float64 value);
  
}


#endif
