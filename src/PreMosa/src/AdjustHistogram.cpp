//
//  AdjustHistogram.cpp
//  CMakeProject
//
//  Created by blasse on 1/18/13.
//
//

#include "AdjustHistogram.h"

#include <iostream>

extern "C" {
#include "image.h"
#include "histogram.h"
}

#include "HelperFunctions.h"
#include "ImageClass.h"
#include "ShowArray.h"
#include "Develop.h"

using namespace Parameter;

namespace AdjustHistogram {
  
  Array * HistogramFitting (Array * projection, ImgParameter & imgPara, Array * reference){
    
    Histogram * h_projection = Make_Histogram(UVAL,256,ValU(1),ValU(0));
    Histagain_Array(h_projection, projection, 0);
    int8 mean_projection = Histogram_Mean(h_projection);
    
    Histogram * h_reference = Make_Histogram(UVAL,256,ValU(1),ValU(0));
    Histagain_Array(h_reference, reference, 0);
    int8 mean_reference = Histogram_Mean(h_reference);
    
    
    
    // Create lookup mapping table
    
    Array * correctedProjection = Copy_Array(projection);
    uint8 * projVals = AUINT8(correctedProjection);
    int * lookup = new int[256];
    
    for (int p=0; p<256; p++) {

      double percent_proj = Bin2Percentile(h_projection, p);      // Match according to the percentiles
      int bin_ref = Percentile2Bin(h_reference, percent_proj);
      lookup[p] = bin_ref;
    
      /*   // Cut brightes/darkest 5% of the values
      if (percent_proj >= 0.95) {
        lookup[p] = 0;
      } else {
        int bin_ref = Percentile2Bin(h_reference, percent_proj);
        lookup[p] = bin_ref;
      }*/
    }
     
    // Adjust each pixel accoring to the lookup table
     
    for (int p=0; p<correctedProjection->size; p++) {
      int projValue = (int) projVals[p];
      projVals[p] = (uint8)lookup[projValue];
    }

    // Cut brightes/darkest 5% of the values
    
    
    //Scale_Array_To_Range(correctedProjection, ValU(0), ValU(255));
    
    
    
    Free_Histogram(h_projection);
    Free_Histogram(h_reference);
    delete [] lookup;
    return correctedProjection;
  }
  
  Array * SmoothLowIntensityRegions (Array * projection, ImgParameter & imgPara){
    
    uint8 * proj    = AUINT8(projection);
    Array  * filteredProjection = Make_Array_With_Shape(PLAIN_KIND, UINT8_TYPE, Coord2(projection->dims[1], projection->dims[0]));
    uint8 * filt    = AUINT8(filteredProjection);
    
    Use_Extend_Boundary();
    
    Frame * f = Make_Frame(projection,Coord2(3,3),Coord2(1, 1));
    Histogram * h = Make_Histogram(UVAL,256,ValU(1),ValU(0));
    
    
    Place_Frame(f,0);
    for (Indx_Type p = 0; p < projection->size; p++){
      
      Empty_Histogram(h);
      Histagain_Array(h,f,0);
      if (proj[p] - Percentile2Bin(h,.5) >= 1 && proj[p] < 20) {
        filt[p] = (uint8)Percentile2Bin(h,.5);
      } else {
        filt[p] = proj[p];
      }
      Move_Frame_Forward(f);
    }
    
    Kill_Histogram(h);
    Kill_Frame(f);
    
    return filteredProjection;
  }
  
  
  
  int Percentile2SmallestBin (Histogram * h, double percentile) {
    
    int currNumber = 0;
    
    for (int i = 0; i < h->nbins; i++) {
      currNumber += h->counts[i];
      
      if (currNumber/h->total >= percentile) {
        return i;
      }
    }

    return h->nbins-1;
  }
  
  
  double Bin2PercentileReversed (Histogram * h, int b) {
    
    int currNumber = 0;
    
    for (int i = 0; i <= b; i++) {
      currNumber += h->counts[i];
    }
    
    return (double)currNumber/h->total;
  }

  
  
  Array * ImproveBorderRegions (Array * image, ImgParameter & imgPara, Array * reference) {
    
    int height = image->dims[1];
    int width = image->dims[0];
    
    
    if (height <= 50 or width <= 50) {
      
      return image;
      
    } else {
      
      uint8 * imgVals = AUINT8(image);
      Array * correctedImg = Copy_Array(image);
      uint8 * correctedVals = AUINT8(correctedImg);
      
      
      // Reference Histogram
      
      Frame * f_reference = Make_Frame(image,Coord2(100,100),Coord2(50, 50));
      Histogram * h_reference = Make_Histogram(UVAL,256,ValU(1),ValU(0));
      Indx_Type center = Coord2IdxA(image, Coord2(ceil(height/2), ceil(width/2)));
      Place_Frame(f_reference, center);
      Histagain_Array(h_reference, f_reference, 0);

      /*
      Histogram * h_reference = Make_Histogram(UVAL,256,ValU(1),ValU(0));
      Histagain_Array(h_reference, reference, 0);
      */
      
      /*
      // Adjust upper border
      
      Frame * f = Make_Frame(image,Coord2(50,width),Coord2(0, 0));
      Histogram * h = Make_Histogram(UVAL,256,ValU(1),ValU(0));
      
      Place_Frame(f, 0);
      Histagain_Array(h, f, 0);
      
      for (Dimn_Type y = 0; y < 50; y++) {
        for (Dimn_Type x = 0; x < width; x++) {
          
          Indx_Type i = Coord2IdxA(image, Coord2(y, x));
          
          double realPercentage = Bin2Percentile(h, imgVals[i]);
          uint8 correctedBin = Percentile2Bin(h_reference, realPercentage);
          
          correctedVals[i] = correctedBin;
        }
      }
      Free_Frame(f);
      Free_Histogram(h);
      
      
      
      // Adjust lower border
      
      f = Make_Frame(image,Coord2(height-50,width),Coord2(height-50, 0));
      h = Make_Histogram(UVAL,256,ValU(1),ValU(0));
      
      Place_Frame(f, 0);
      Histagain_Array(h, f, 0);
      
      for (Dimn_Type y = height-50; y < height; y++) {
        for (Dimn_Type x = 0; x < width; x++) {
          
          Indx_Type i = Coord2IdxA(image, Coord2(y, x));
          
          double realPercentage = Bin2Percentile(h, imgVals[i]);
          uint8 correctedBin = Percentile2Bin(h_reference, realPercentage);
          correctedVals[i] = correctedBin;
        }
      }
      Free_Frame(f);
      Free_Histogram(h);
      
      // Adjust left border
      
      Frame *f = Make_Frame(image,Coord2(height-100,50),Coord2(50, 0));
      Histogram *h = Make_Histogram(UVAL,256,ValU(1),ValU(0));
      
      Place_Frame(f, 0);
      Histagain_Array(h, f, 0);
      
      for (Dimn_Type y = 50; y < height-50; y++) {
        for (Dimn_Type x = 0; x < 50; x++) {
          
          Indx_Type i = Coord2IdxA(image, Coord2(y, x));
          
          double realPercentage = Bin2Percentile(h, imgVals[i]);
          uint8 correctedBin = Percentile2Bin(h_reference, realPercentage);
          correctedVals[i] = correctedBin;
        }
      }
      Free_Frame(f);
      Free_Histogram(h);

      */
      // Adjust right border
      
      Frame *f = Make_Frame(image,Coord2(height-100,50),Coord2(50, width-50));
      Histogram *h = Make_Histogram(UVAL,256,ValU(1),ValU(0));
      
      Place_Frame(f, 0);
      Histagain_Array(h, f, 0);
      
      for (Dimn_Type y = 50; y < height-50; y++) {
        for (Dimn_Type x = width-50; x < width; x++) {
          
          Indx_Type i = Coord2IdxA(image, Coord2(y, x));
          
          double realPercentage = Bin2PercentileReversed(h, (int)imgVals[i]);
          correctedVals[i] = (uint8) Percentile2SmallestBin(h_reference, realPercentage);
          
          std::cout << realPercentage << "   " << (int) imgVals[i] << " - " << Percentile2SmallestBin(h_reference, realPercentage) << std::endl;
        }
      }
      Free_Frame(f);
      Free_Histogram(h);

    
      
      Write_Image("Corrected.tif", correctedImg, DONT_PRESS);
      
      return correctedImg;
    }
        
  }
  
  
}