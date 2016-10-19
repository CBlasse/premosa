//
//  LocalIntegrals.cpp
//  CMakeProject
//
//  Created by blasse on 1/11/13.
//
//

#include "LocalIntegrals.h"

using namespace HelperFunctions;


namespace Locally {
  
  
  Array * GetSDImg (Array * image) {
    
    Array  * sdImg = Make_Array_With_Shape(PLAIN_KIND,UINT32_TYPE,Coord2(image->dims[1], image->dims[0]));
    uint32 * sdVals    = AUINT32(sdImg);
    
    Use_Reflective_Boundary();
    /*
     #pragma omp parallel for
     for (Indx_Type p = 0; p < image->size; ++p)
     {
     Frame* f = Make_Frame (image, Coord2(60, 60), Coord2 (30, 30));
     Histogram* h = Make_Histogram (UVAL, 256, ValU (1), ValU (0));
     
     Place_Frame (f, p);
     
     Histagain_Array (h, f, 0);
     
     sdVals[p] = (uint32) Histogram_Sigma (h);
     
     Free_Frame (f);
     Free_Histogram (h);
     }
     */
    
    Use_Reflective_Boundary();
    
    Frame * f = Make_Frame(image,Coord2(40,40),Coord2(20,20));
    Histogram * h = Make_Histogram(UVAL, 256,ValU(1),ValU(0));
    Place_Frame(f,0);

    
    Indx_Type lp = 0;
    for (Dimn_Type y = 0; y < image->dims[1]; y += 10){
      for (Dimn_Type x = 0; x < image->dims[0]; x += 10){
        Place_Frame(f, Coord2IdxA(image, Coord2(y, x)));
        
        Empty_Histogram(h);
        Histagain_Array(h,f,0);
        
        uint32 sd = (uint32) Histogram_Sigma(h);
        
        for (Dimn_Type d_1 = -5; d_1 <= 5 ; d_1++) {
          for (Dimn_Type d_2 = -5; d_2 <= 5 ; d_2++) {

            if (y+d_1 >=0 and x+d_2 >= 0 and y+d_1 < image->dims[1] and x+d_2 < image->dims[0]) {
              Indx_Type p = Coord2IdxA(image, Coord2(y+d_1, x+d_2));
              sdVals[p] = sd;
            }
          }
        }
        Move_Frame_Forward(f);
      }
    }
    
    /*
    
    
    Histogram * h = Make_Histogram(UVAL, 256,ValU(1),ValU(0));
    Place_Frame(f,0);
    
    for (Indx_Type p=0; p<image->size; p=p+1) {
      Empty_Histogram(h);
      Histagain_Array(h,f,0);
      
      sdVals[p] = (uint32) Histogram_Sigma(h);
      
      Coordinate * co = Idx2CoordA(image, p);
      
      Move_Frame_Forward(f);
      
    }*/
    
    Write_Image("Sd.tif", sdImg, DONT_PRESS);
    
    Kill_Frame(f);
    Kill_Histogram(h);
    
    return sdImg;
  }
  
  
  Array * GetAvgImg (Array * image) {
    
    Array  * avgImg = Make_Array_With_Shape(PLAIN_KIND,UINT32_TYPE,Coord2(image->dims[1], image->dims[0]));
    uint32 * avgVals    = AUINT32(avgImg);
    
    Use_Reflective_Boundary();
        
    Frame * f = Make_Frame(image,Coord2(40,40),Coord2(20,20));
    Histogram * h = Make_Histogram(UVAL, 256,ValU(1),ValU(0));
    Place_Frame(f,0);
    
    
    Indx_Type lp = 0;
    for (Dimn_Type y = 0; y < image->dims[1]; y += 10){
      for (Dimn_Type x = 0; x < image->dims[0]; x += 10){
        Place_Frame(f, Coord2IdxA(image, Coord2(y, x)));
        
        Empty_Histogram(h);
        Histagain_Array(h,f,0);
        
        uint32 avg = (uint32) Histogram_Mean(h);
        
        for (Dimn_Type d_1 = -5; d_1 <= 5 ; d_1++) {
          for (Dimn_Type d_2 = -5; d_2 <= 5 ; d_2++) {
            
            if (y+d_1 >=0 and x+d_2 >= 0 and y+d_1 < image->dims[1] and x+d_2 < image->dims[0]) {
              Indx_Type p = Coord2IdxA(image, Coord2(y+d_1, x+d_2));
              avgVals[p] = avg;
            }
          }
        }
        Move_Frame_Forward(f);
      }
    }
    
    Write_Image("Avg.tif", avgImg, DONT_PRESS);   
    
    Kill_Frame(f);
    Kill_Histogram(h);
    
    return avgImg;
  }

  
  Array * LocalMedianFilter (Array * image){
    
    Array  * filteredImg = Copy_Array(image);
    uint8 * filteredVals    = AUINT8(filteredImg);

    Array * img2 = Convert_Array(image, RGB_KIND, UINT8_TYPE, 8, 0);
    static Color_Bundle MAGENTA = { SET_PIX, {255}, {0}, {255}, {0} };
    
    
    Use_Reflective_Boundary();
    
    Frame * f = Make_Frame(image,Coord2(3,3),Coord2(1,1));
    Histogram * h = Make_Histogram(UVAL, 256,ValU(1),ValU(0));
    Place_Frame(f,0);
    
    Indx_Type lp = 0;
    for (Dimn_Type y = 0; y < image->dims[1]; y++){
      for (Dimn_Type x = 0; x < image->dims[0]; x++){
        Place_Frame(f, Coord2IdxA(image, Coord2(y, x)));
        
        Empty_Histogram(h);
        Histagain_Array(h,f,0);
        
        uint32 sd = (uint32) Histogram_Sigma(h);
        uint8 median = (uint8) Percentile2Bin(h,.5);
        
        if (sd < 3.5 && median <= 80) {
          Indx_Type p = Coord2IdxA(image, Coord2(y, x));
          filteredVals[p] = median;
          Draw_Point(img2, &MAGENTA, Coord2(y, x));
        }
        Move_Frame_Forward(f);
      }
    }
    Write_Image("MarkedMedianFilteredImg.tif", img2, DONT_PRESS);
    Write_Image("MedianFilteredImg.tif", filteredImg, DONT_PRESS);
    
    Kill_Frame(f);
    Kill_Histogram(h);
    
    Free_Array(img2);
    
    
    return filteredImg;

    
  }
  
}
