////
////  FilterFFTW.cpp
////  CMakeProject
////
////  Created by blasse on 1/28/13.
////
////
//
//#include "FilterFFTW.h"
//
//
//#include <iostream>
//#include <ostream>
//#include <cmath>
//#define _USE_MATH_DEFINES
////#include "fftw3.h"
//
//
//extern "C" {
//#include "image.h"
//}
//#include "ImageClass.h"
//#include "Parameter.h"
//#include "HelperFunctions.h"
//#include "Develop.h"
//
//
//using namespace HelperFunctions;
//
//
//namespace FilteringFFTW {
//  
//  int Next_Power_Of_2(int m) { 
//    int n = 2;
//    while (n < m)
//      n <<= 1;
//    return (n);
//  }
//  
//  
//  Array * FFTPreprocessing (ImageClass::Image & img) {
//    
//    int maxLength = fmax(img.GetHeight(), img.GetWidth());
//    int newLength = Next_Power_Of_2(3*maxLength);
//    
//    Use_Reflective_Boundary();
//    Array * paddedImage = Pad_Array(img.GetImage(), Coord2(ceil((newLength-img.GetHeight())/2), ceil((newLength-img.GetWidth())/2)), Coord2(newLength, newLength));
//    uint8 * paddedVals = AUINT8(paddedImage);
//    Array * output = Make_Array_With_Shape(PLAIN_KIND, UINT8_TYPE, Coord2(newLength, newLength));
//    uint8 * outputVals = AUINT8(output);
//    
//    int threshold = sqrt(maxLength/2 * maxLength/2 + maxLength/2 * maxLength/2);
//    int maxThreshold = sqrt(newLength/2 * newLength/2 + newLength/2 * newLength/2);
//    
//    double scale = newLength/(newLength/maxLength);
//    
//    for (Dimn_Type x = 0; x < newLength; x++) {
//      for (Dimn_Type y = 0; y < newLength; y++) {
//        
//        float xDist = abs((ceil(newLength/2)) - x);
//        float yDist = abs((ceil(newLength/2)) - y);
//        float distance = sqrt(xDist*xDist + yDist*yDist);
//        
//        Indx_Type p = Coord2IdxA(output, Coord2(y, x));
//        outputVals[p] = (uint8) (exp(-(xDist*xDist + yDist*yDist)/(2*scale*scale)) * paddedVals[p]);
//      }
//    }
//    
//    
//    Write_Image("Prepared.tif", output, DONT_PRESS);
//    //Write_Image("Padded2.tif", paddedImage, DONT_PRESS);
//
//    return output;
//  }
//  
//  
//  Array * FFTFilter_Quadrupled(ImageClass::Image & img, Parameter::ImgParameter & imgPara){
//    
//    Array * paddedImage = FFTPreprocessing (img);
//    int newLength = paddedImage->dims[0];
//    int maxOldLength = fmax((double) img.GetHeight(), (double) img.GetWidth());
//    
//    uint8 * imgVals = AUINT8(paddedImage);
//    int imgSize = newLength * newLength;
//    
//    
//    double cutoff = ceil(newLength/3);
//    
//    
//    fftw_complex * data         = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * imgSize);
//    fftw_complex * fft_result   = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * imgSize);
//    fftw_complex * invertedData = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * imgSize);
//    
//    fftw_plan plan_forward = fftw_plan_dft_1d(imgSize, data, fft_result, FFTW_FORWARD, FFTW_ESTIMATE);
//    fftw_plan plan_backward = fftw_plan_dft_1d(imgSize, fft_result, invertedData, FFTW_BACKWARD, FFTW_ESTIMATE);
//    
//#ifdef DEVELOP
//    Array * fftImage = Make_Array_With_Shape(PLAIN_KIND, FLOAT32_TYPE, Coord2(newLength, newLength));
//    float32 * fftVals = AFLOAT32(fftImage);
//#endif
//    
//    for (int i=0; i<imgSize; i++) {         // Assignment of the data values
//      data[i][0] = imgVals[i];
//      data[1][1] = 0.0;
//    }
//    
//    fftw_execute(plan_forward);             // FFT
//    
//    for (Dimn_Type y=0; y<newLength; y++) {
//      for (Dimn_Type x=0; x<newLength; x++) {
//        
//        Indx_Type i = Coord2IdxA(paddedImage, Coord2(y, x));
//        Indx_Type iSwitched;
//        int xSwitched, ySwitched;
//        
//        //std::cout << fft_result[i][0]/x << " - " << fft_result[i][1]/y << std::endl;
//        
//        if (x > 100) {        // TODO: Why does this shift occur?
//          xSwitched = (x+(newLength/2))%newLength;
//          ySwitched = (1+y+(newLength/2))%newLength;
//          iSwitched = Coord2IdxA(paddedImage, Coord2(ySwitched, xSwitched));
//          
//        } else {
//          xSwitched = (x+(newLength/2))%newLength;
//          ySwitched = (y+(newLength/2))%newLength;
//          iSwitched = Coord2IdxA(paddedImage, Coord2(ySwitched, xSwitched));
//        }
//        
//        float amplitude = sqrt(fft_result[i][0]*fft_result[i][0] + fft_result[i][1]*fft_result[i][1]);
//        float val;
//        
//        float xDist = abs((ceil(newLength/2)) - xSwitched);
//        float yDist = abs((ceil(newLength/2)) - ySwitched);
//        float distance = sqrt(xDist*xDist + yDist*yDist);
//        
//        if (distance < cutoff) {
//          fft_result[i][0] = fft_result[i][0];
//          fft_result[i][1] = fft_result[i][1];
//          //fftVals[iSwitched] = log(amplitude);
//          
//          //fft_result[i][0] = fft_result[i][0]*exp(-0.5*(distance/cutoff));
//          //fft_result[i][1] = fft_result[i][1]*exp(-0.5*(distance/cutoff));
//#ifdef DEVELOP          
//          fftVals[iSwitched] = log(amplitude* exp(-0.5*(distance/cutoff)));
//#endif
//          
//        } else {
//          
//          fft_result[i][0] = 0;
//          fft_result[i][1] = 0;
//#ifdef DEVELOP
//          fftVals[iSwitched] = 0;
//#endif
//        }
//        
//        
//      }
//      
//      
//    }
//    
//#ifdef DEVELOP
//    Scale_Array_To_Range(fftImage, ValF32(0.0), ValF32(255.0));
//    
//    Array * convertedArray = Convert_Array(paddedImage, PLAIN_KIND, UINT8_TYPE, 8, 0);
//    uint8 * convVals = AUINT8(convertedArray);
//    
//    for (int i=0; i<imgSize; i++) {
//      convVals[i] = (uint8) fftVals[i];
//      //convVals[i] = (uint8) log(fftVals[i]/imgSize);
//    }
//    
//    
//    Write_Image("FFT.tif", convertedArray, DONT_PRESS);
//#endif
//    
//    fftw_execute(plan_backward);
//    
//    
//    Array * imgForScaling = Make_Array_With_Shape(PLAIN_KIND, INT32_TYPE, Coord2(newLength, newLength));
//    int32 * valForScaling = AINT32(imgForScaling);
//    
//    for(int i = 0 ; i < imgSize ; i++ ) {
//      valForScaling[i] = invertedData[i][0]/imgSize;
//    }
//    
//    
//    Range_Bundle range;
//    Array_Range(&range, imgForScaling);
//    
//    std::cout << range.minval.ival << " - " << range.maxval.ival << std::endl;
//    
//    //Scale_Array_To_Range(imgForScaling, ValI(10), ValI(200));
//    Scale_Array_To_Range(imgForScaling, ValI(0), ValI(255));
//    
//    
//    Array * invertedImage = Copy_Array(paddedImage);
//    uint8 * invImgVals = AUINT8(invertedImage);
//    
//    for(int i = 0 ; i < imgSize ; i++ ) {
//      /*
//       if ((uint8)valForScaling[i] <= 0) {
//       invImgVals[i] = 0;
//       } else if ((uint8)valForScaling[i] >= 255){
//       invImgVals[i] = 255;
//       } else {
//       invImgVals[i] = (uint8)valForScaling[i];
//       }*/
//      invImgVals[i] = (uint8)valForScaling[i];
//    }
//    
//    Clip_Array_Inplace(invertedImage, Coord2(ceil((newLength-img.GetHeight())/2), ceil((newLength-img.GetWidth())/2)), Coord2(floor((newLength+img.GetHeight())/2)-1, floor((newLength+img.GetWidth())/2)-1));
//    Scale_Array_To_Range(invertedImage, ValU(0), ValU(255));
//    
//    Write_Image("invertedData.tif", invertedImage, DONT_PRESS);
//    
//    fftw_free(data);
//    fftw_free(fft_result);
//    fftw_free(invertedData);
//    fftw_destroy_plan(plan_forward);
//    fftw_destroy_plan(plan_backward);
//    
//    Free_Array(imgForScaling);
//    Free_Array(paddedImage);
//    
//    return invertedImage;
//    
//  }
//}