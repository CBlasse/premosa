//
//  VectorFunctions.h
//
//  Created by Corinna Blasse on 04/01/15.
//
//

#ifndef VECTORFUNCTIONS_H
#define VECTORFUNCTIONS_H

#include <vector>

namespace VectorFunctions {
  
  template <typename T> void PrintVector (std::vector<T> &v);
  
  template <typename T> double VectorMean (std::vector<T> &v);
  template <typename T> double VectorMedian (std::vector<T> &v);
  template <typename T> double VectorVariance (std::vector<T> &v);
  template <typename T> double VectorSd (std::vector<T> &v);
  template <typename T> double VectorEntropy (std::vector<T>& v);
}


#endif
