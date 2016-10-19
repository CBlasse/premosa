//
//  FilterFFTW.h
//  CMakeProject
//
//  Created by blasse on 1/28/13.
//
//

#ifndef FILTERFFTW_H
#define FILTERFFTW_H


namespace Parameter {
  struct ImgParameter;
}
namespace ImageClass {
  class Image;
}

//struct Array;

extern "C"{
#include "array.h"
}




namespace FilteringFFTW {
  
  int Next_Power_Of_2(int m);

  Array * FFTPreprocessing (ImageClass::Image & img);
  
  Array * FFTFilter_Quadrupled(ImageClass::Image & img, Parameter::ImgParameter & imgPara);
  
}

#endif 