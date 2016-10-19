//
//  ZProjection.cpp
//
//  Created by Corinna Blasse on 25.07.12.
//

#include <iostream>

#include "ZProjection.h"

extern "C" {
#include "histogram.h"
}

#include "MylibValues.h"

using namespace MylibValues;


namespace ProjectionMethod {
  
  /***************
   * Computes the max z-projection of an image stack (8-bit)
   ***************/
  
  Array * MaxProjection_8bit (RasterizedImage & img, ParaSet & imgPara){
    Indx_Type q;
    int       zmax;
    
    Array  *maxim = Make_Array(PLAIN_KIND,UINT8_TYPE,2, img.GetImage()->dims);
    uint8 * max   = AUINT8(maxim);
    uint8 * val   = AUINT8(img.GetImage());
    
    for (Indx_Type p = 0; p < img.GetArea(); p++){ 
      zmax = val[p];
      q = p + img.GetArea();
      
      for (Dimn_Type d = 1; d < img.GetDepth(); d++){ 
        if (val[q] >= zmax) {
          zmax = val[q];

        }
        q += img.GetArea();
      }
      max[p] = zmax;
    }
    
#ifdef DEVELOP
    ShowArray(img, imgPara, "maxim.tif",maxim);
#endif
 
    if (imgPara.verbose){ 
      std::cout << "Max-projection computed" << std::endl;
    }
        
    return (maxim);
  }
  
  /***************
   * Computes the max z-projection of an image stack (16-bit)
   ***************/
  
  Array * MaxProjection_16bit (RasterizedImage & img, ParaSet & imgPara){
    Indx_Type q;
    int       zmax;
    
    Array  *maxim = Make_Array(PLAIN_KIND,UINT16_TYPE,2, img.GetImage()->dims);
    uint16 * max   = AUINT16(maxim);
    uint16 * val   = AUINT16(img.GetImage());
    
    for (Indx_Type p = 0; p < img.GetArea(); p++){
      zmax = val[p];
      q = p + img.GetArea();
      
      for (Dimn_Type d = 1; d < img.GetDepth(); d++){
        if (val[q] >= zmax) {
          zmax = val[q];
          
        }
        q += img.GetArea();
      }
      max[p] = zmax;
    }
    
#ifdef DEVELOP
    ShowArray(img, imgPara, "maxim.tif",maxim);
#endif
    
    if (imgPara.verbose){
      std::cout << "Max-projection computed" << std::endl;
    }
    
    return (maxim);
  }

  
  
  /***************
  *
  * Computes the max z-projection of an image using a subset of the image stack 
  *   defined by a heiht map and a distance value
  * 8-bit Images
  *
  ***************/
  
  Array * Substack_UsingMask_8bit(RasterizedImage & img, ParaSet & imgPara, Array * median, bool onlyTopPixel){
    
    Array * substack = Make_Array(PLAIN_KIND, UINT8_TYPE, 3, img.GetImage()->dims);
    Array_Op_Scalar(substack, SET_OP, UINT8_TYPE, ValU(0));
      
    uint8   *maxSt = AUINT8(substack);
    uint8   *val = AUINT8(img.GetImage());
    uint32  *map = AUINT32(median); 
    
    
    Indx_Type p = 0;
    
    int stackSize = std::min(2*imgPara.distance+1, img.GetDepth());
    
    Use_Extend_Boundary();
    Frame * f = Make_Frame(img.GetImage(), Coord3(stackSize, imgPara.radius, imgPara.radius), Coord3(ceil((double)stackSize/(double)2), 0, 0));
    Histogram * h = Make_Histogram(UVAL,256,ValU(1),ValU(0));
    Place_Frame(f, 0);
    
    
    uint32      z, zmax;
    Size_Type   *hcnt;
    int         sum, thr;
    Indx_Type   tmpP;
    uint8       maxVal;

    
    for (Dimn_Type y=0; y < img.GetGridY(); y++) {
      for (Dimn_Type x=0; x < img.GetGridX(); x++) {
        
        p = y*img.GetGridX() + x;
        z = map[p];
        
        p = Coord2IdxA(img.GetImage(), Coord3(z, y*imgPara.radius, x*imgPara.radius));
        Place_Frame(f, p);
        Histagain_Array(Empty_Histogram(h),f,0);
        hcnt = h->counts;
        
        sum = 0;
        for (thr = 255; thr > 0; thr--)
        { 
          sum += hcnt[thr];
          if (sum >= imgPara.threshold*img.GetDepth())
            break;
        }
        
        for (int r=0; r < imgPara.radius; r++) {
          for (int s=0; s < imgPara.radius; s++) {
            
            Indx_Type yCoord = (Indx_Type) y*imgPara.radius + r;
            Indx_Type xCoord = (Indx_Type) x*imgPara.radius + s;
              
            p = Coord2IdxA(img.GetImage(), Coord3(z, yCoord, xCoord));
            maxVal = val[p];
            zmax = z;
            
            for (int d = -imgPara.distance; d <= imgPara.distance; d++) {
              tmpP = p + d*img.GetArea();
              
              if ((tmpP > imgPara.layer) && (tmpP < img.GetDepth()*img.GetArea())) {
                if (onlyTopPixel == true) {
                  if ((val[tmpP] > maxVal) && (val[tmpP] >= thr)) {
                    maxVal = val[tmpP];
                    zmax = z + d;
                    
                  }
                } else {
                  p = Coord2IdxA(substack, Coord3(z+d, yCoord, xCoord));
                  maxSt[p] = (uint8)val[tmpP];
                }
                
              }            
            }
            if (onlyTopPixel == true) {
              p = Coord2IdxA(substack, Coord3(zmax, yCoord, xCoord));
              maxSt[p] = (uint8)maxVal;
            }
            
          }
        }
      }
    }
    
    Free_Frame(f);
    Free_Histogram(h);
    
    
#ifdef DEVELOP
    ShowArray(img, imgPara, "substacks.tif",substack);
#endif
    
    if (imgPara.verbose){ 
      std::cout << "Restricted max-projection computed" << std::endl; 
    }
    
    return substack;
  }
  
  
  /***************
   *
   * Computes the max z-projection of an image using a subset of the image stack
   *   defined by a heiht map and a distance value
   *  16-bit Images
   *
   ***************/
  
  Array * Substack_UsingMask_16bit(RasterizedImage & img, ParaSet & imgPara, Array * median, bool onlyTopPixel){
    
    uint32   *map = AUINT32(median);
    
    /******
     * Speed up by using an 8-bit image to determine the layer with the highest variance
     ******/
    Array * substack = Make_Array(PLAIN_KIND, UINT8_TYPE, 3, img.GetImage()->dims);
    Array_Op_Scalar(substack, SET_OP, UINT8_TYPE, ValU(0));
    uint8   *maxSt = AUINT8(substack);

    Array * img_8bit = Copy_Array(img.GetImage());      //something is wrong here
    Scale_Array_To_Range(img_8bit,ValU(0),ValU(255));
    Convert_Array_Inplace(img_8bit, PLAIN_KIND,UINT8_TYPE,8,0);
    uint8   *val = AUINT8(img_8bit);
    
    /******
     * Otherwise
     ******/
    /*
    Array * substack = Make_Array(PLAIN_KIND, UINT16_TYPE, 3, img.GetImage()->dims);
    uint16   *maxSt = AUINT16(substack);
    uint16   *val = AUINT16(img.GetImage());
    */
    
    Indx_Type p;
        
    int stackSize = std::min(2*imgPara.distance+1, img.GetDepth());
    
    Use_Extend_Boundary();
    Frame * f = Make_Frame(img.GetImage(), Coord3(stackSize, imgPara.radius, imgPara.radius), Coord3(ceil((double)stackSize/(double)2), 0, 0));
    Histogram * h = Make_Histogram(UVAL,256,ValU(1),ValU(0));
    Place_Frame(f, 0);
    
    
    uint32      z, zmax;
    Size_Type   *hcnt;
    int         sum, thr;
    Indx_Type   tmpP;
    uint8       maxVal;
    
    Range_Bundle range;
    Array_Range(&range, median);
    std::cout << range.minval.uval << " " << range.maxval.uval << std::endl;
    
    //std::cout << map[(Indx_Type)1466] << std::endl;
    
    for (Dimn_Type y=0; y < img.GetGridY(); y++) {
      for (Dimn_Type x=0; x < img.GetGridX(); x++) {
        
        //p = Coord2IdxA(median, Coord2(y, x));
        p = (Indx_Type)y*img.GetGridX() + x;
        
        z = map[p];     // TODO: Cheating
       
        std::cout << "P " << p << "    " << map[1466] << " " << x << ", " << y << " , " << z <<  std::endl;
        p = Coord2IdxA(img.GetImage(), Coord3(z, y*imgPara.radius, x*imgPara.radius));
        Place_Frame(f, p);
        Histagain_Array(Empty_Histogram(h),f,0);
        hcnt = h->counts;
        
        sum = 0;
        for (thr = 255; thr > 0; thr--)
        {
          sum += hcnt[thr];
          if (sum >= imgPara.threshold*img.GetDepth())
            break;
        }
        
        for (int r=0; r < imgPara.radius; r++) {
          for (int s=0; s < imgPara.radius; s++) {
            
            p = Coord2IdxA(img.GetImage(), Coord3(z, y*imgPara.radius+r, x*imgPara.radius+s));
            maxVal = val[p];
            zmax = z;
            
            for (int d = -imgPara.distance; d <= imgPara.distance; d++) {
              tmpP = p + d*img.GetArea();
              
              if ((tmpP > 0) && (tmpP < img.GetDepth()*img.GetArea())) {
                if (onlyTopPixel == true) {
                  if ((val[tmpP] > maxVal) && (val[tmpP] >= thr)) {
                    maxVal = val[tmpP];
                    zmax = z + d;
                    
                  }
                } else {
                  p = Coord2IdxA(substack, Coord3(z+d, y*imgPara.radius+r, x*imgPara.radius+s));
                  maxSt[p] = val[tmpP];
                }
                
              }
            }
            if (onlyTopPixel == true) {
              p = Coord2IdxA(substack, Coord3(zmax, y*imgPara.radius+r, x*imgPara.radius+s));
              maxSt[p] = maxVal;
            }
            
          }
        }
      }
    }
    
    Free_Frame(f);
    Free_Histogram(h);
    
    
#ifdef DEVELOP
    ShowArray(img, imgPara, "substacks.tif",substack);
#endif
    
    if (imgPara.verbose){
      std::cout << "Restricted max-projection computed" << std::endl;
    }
    
    return substack;
  }

  
  
  
  /***************
   * Compute the max z-projection of image in maxim, and record the z-depth of 
   *   each max value in index.
   ***************/
  
  Array *MaxIndices_8bit(RasterizedImage & img, ParaSet & imgPara, Array *maxim){
    
    Indx_Type q;
    int       best, zmax;
    
    Array  *index = Make_Array(PLAIN_KIND,UINT32_TYPE,2,img.GetImage()->dims);
    uint32 *idx   = AUINT32(index);
    uint8  *val   = AUINT8(img.GetImage());
    uint8  *max   = AUINT8(maxim);
    
    for (Indx_Type p = 0; p < img.GetArea(); p++){ 
      zmax = max[p];
      q = p;
      best = 0;
      for (Dimn_Type d = 0; d < img.GetDepth(); d++){ 
        if (val[q] == zmax){
          best = d;
        }
        q += img.GetArea();
      }
      idx[p] = best;
    }
#ifdef DEVELOP
    ShowArray(img, imgPara, "index.tif",index);
#endif
    
    if (imgPara.verbose){ 
      std::cout << "Indices of max projection computed" << std::endl; 
    }
    
    return (index);
  }
  
  
  /***************
   * Compute the max z-projection of image in maxim, and record the z-depth of
   *   each max value in index.
   ***************/
  
  Array *MaxIndices_16bit(RasterizedImage & img, ParaSet & imgPara, Array *maxim){
    
    Indx_Type q;
    int       best, zmax;
    
    Array  *index = Make_Array(PLAIN_KIND,UINT32_TYPE,2,img.GetImage()->dims);
    uint32 *idx   = AUINT32(index);
    uint16  *val   = AUINT16(img.GetImage());
    uint16  *max   = AUINT16(maxim);
    
    for (Indx_Type p = 0; p < img.GetArea(); p++){
      zmax = max[p];
      q = p;
      for (Dimn_Type d = 0; d < img.GetDepth(); d++){
        if (val[q] == zmax){
          best = d;
        }
        q += img.GetArea();
      }
      idx[p] = best;
    }
    
#ifdef DEVELOP
    ShowArray(img, imgPara, "index.tif",index);
#endif
    
    if (imgPara.verbose){
      std::cout << "Indices of max projection computed" << std::endl;
    }
    
    return (index);
  }


  Array * GetSDImg (Array * image, int radius) {
    
    Array  * sdImg = Make_Array_With_Shape(PLAIN_KIND,UINT32_TYPE,Coord2(image->dims[1], image->dims[0]));
    uint32 * sdVals    = AUINT32(sdImg);
    
    double xPad = ceil(image->dims[0]/radius)* radius + 1;
    double yPad = ceil(image->dims[1]/radius)* radius + 1;
    
    if (radius % 2 > 0) {
      radius++;
    }
    
    Array * paddedImg = Pad_Array(image, Coord2(0, 0), Coord2(yPad, xPad));
    
    Use_Reflective_Boundary();
    
    Frame * f = Make_Frame(paddedImg,Coord2(radius,radius),Coord2(0,0));
    Histogram * h;
    if (image->type == UINT8_TYPE) {
      h = Make_Histogram(UVAL, 256,ValU(1),ValU(0));
    } else {
      h = Make_Histogram(UVAL, 65536,ValU(1),ValU(0));
    }
    Place_Frame(f,0);
    
    for (Dimn_Type y = 0; y < paddedImg->dims[1]; y += radius){
      for (Dimn_Type x = 0; x < paddedImg->dims[0]; x += radius){
        Place_Frame(f, Coord2IdxA(paddedImg, Coord2(y, x)));
        
        Empty_Histogram(h);
        Histagain_Array(h,f,0);
        
        uint32 sd = (uint32) Histogram_Sigma(h);
        
        for (Dimn_Type d_1 = -radius/2; d_1 <= radius/2 ; d_1++) {
          for (Dimn_Type d_2 = -radius/2; d_2 <= radius/2 ; d_2++) {
            
            if (y+d_1 >=0 and x+d_2 >= 0 and y+d_1 < image->dims[1] and x+d_2 < image->dims[0]) {
              Indx_Type p = Coord2IdxA(image, Coord2(y+d_1, x+d_2));
              sdVals[p] = sd;
            }
          }
        }
        Move_Frame_Forward(f);
      }
    }
    
    
    Free_Array(paddedImg);
    Kill_Frame(f);
    Kill_Histogram(h);
    
    return sdImg;
  }


}