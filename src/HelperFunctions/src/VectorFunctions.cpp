//
//  VectorFunctions.cpp
//
//  Created by Corinna Blasse on 04/01/15.
//


#include "VectorFunctions.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <unordered_map>

namespace VectorFunctions {
  

  ////////////////////////////////////////////////////////////////////////////////
  /**
   Prints the vector
   
   @param vector
   //////////////////////////*/

  template <typename T>
  void PrintVector (std::vector<T> &v) {
    for (auto & item: v ) {
      std::cout << item << ", ";
    }
    std::cout << std::endl;
  }

  // Instantiate functions for the supported template type parameters
  template void PrintVector<int> (std::vector<int> &v);
  template void PrintVector<double> (std::vector<double> &v);


  
  ////////////////////////////////////////////////////////////////////////////////
  /**
   Computes the average value
   
   @param vector
   //////////////////////////*/

  template <typename T>
  double VectorMean (std::vector<T> &v) {
    
    if (v.size() == 0) {
      return 0;
    }
    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
  }

  // Instantiate functions for the supported template type parameters
  template double VectorMean<int> (std::vector<int> &v);
  template double VectorMean<double> (std::vector<double> &v);



  ////////////////////////////////////////////////////////////////////////////////
  /**
   Computes the median value
   
   @param vector
   //////////////////////////*/

  template <typename T>
  double VectorMedian (std::vector<T> &v) {
    
    if (v.size() == 0) {
      return 0;
    }
    sort(v.begin(), v.end());
    if (v.size() % 2 == 0) {
      return (v[(v.size()-1)/2] + v[(v.size()+1)/2]) /2;
    } else {
      return v[v.size()/2];
    }
  }

  // Instantiate functions for the supported template type parameters
  template double VectorMedian<int> (std::vector<int> &v);
  template double VectorMedian<double> (std::vector<double> &v);

  
  
  ////////////////////////////////////////////////////////////////////////////////
  /**
   Computes the variance of the values
   
   @param vector
   //////////////////////////*/

  template <typename T>
  double VectorVariance (std::vector<T> &v) {
    double mean = VectorMean<T>(v);
    double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
    
    return (sq_sum / v.size() - mean * mean);
  }

  // Instantiate functions for the supported template type parameters
  template double VectorVariance<int> (std::vector<int> &v);
  template double VectorVariance<double> (std::vector<double> &v);

  
  
  ////////////////////////////////////////////////////////////////////////////////
  /**
   Computes the standard deviation of the values
   
   @param vector
   //////////////////////////*/

  template <typename T>
  double VectorSd (std::vector<T> &v) {
    return sqrt(VectorVariance(v));
  }

  // Instantiate functions for the supported template type parameters
  template double VectorSd<int> (std::vector<int> &v);
  template double VectorSd<double> (std::vector<double> &v);


  
  ////////////////////////////////////////////////////////////////////////////////
  /**
   Computes the entropy of the values
   
   @param vector
   //////////////////////////*/

  template <typename T>
  double VectorEntropy (std::vector<T>& v){
    
    std::unordered_map<int, int> freq;
    
    for (const auto& intensity : v) {
      ++freq[intensity];
    }
    
    double entropy = 0;
    
    for (const auto& element : freq) {
      double p = (double)element.second / (double)v.size();
      if (p > 1.e-20) {
        entropy -= (double)p*log(p)/log(2);
      }
    }
    return entropy;
  }

  // Instantiate functions for the supported template type parameters
  template double VectorEntropy<int> (std::vector<int> &v);
  template double VectorEntropy<double> (std::vector<double> &v);


}