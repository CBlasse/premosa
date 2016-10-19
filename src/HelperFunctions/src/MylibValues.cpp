//
//  MylibValues.cpp
//
//  Created by Corinna Blasse on 04/01/15.
//


#include "MylibValues.h"

namespace MylibValues {
  
  ////////////////////////////////////////////////////////////////////////////////
  /**
   Converts a number into a Mylib value
   
   @param value
   //////////////////////////*/

  
  Value ValU(uint64 value){
    Value v;
    v.uval = value;
    return v;
  }
  
  Value ValU8(uint8 value){
    Value v;
    v.uval = value;
    return v;
  }
  
  Value ValI(int64 value){
    Value v;
    v.ival = value;
    return v;
  }
  
  Value ValF(float64 value){
    Value v;
    v.fval = value;
    return v;
  }
  
}