////
////  ZPlaneSelection.cpp
////  proj_temp_cpp
////
////  Created by Corinna Blasse on 26.07.12.
////  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
////
//
//#include <iostream>
//
//#include "ZPlaneSelection.h"
//
//using namespace HelperFunctions;
//using namespace ImageClass;
//using namespace Parameter;
//
//namespace ZPlane {
//  
//  /***********************
//   * Select the z-plane with the most occurrences of bright pixel
//   **********************/
//  
//  Array * PlaneSelection (Image & img, ParaSet & imgPara, Array *maxim, Array *index){
//    
//    int tmpGridX = ceil ((double)img.GetWidth()/(double)imgPara.radius);
//    int tmpGridY  = ceil ((double)img.GetHeight()/(double)imgPara.radius);
//    
//    Array  *level = Make_Array_With_Shape(PLAIN_KIND,UINT32_TYPE,Coord2(tmpGridY,tmpGridX));
//    uint32 *lev   = AUINT32(level);
//    
//    uint32 *idx   = AUINT32(index);
//    uint8  *max   = AUINT8(maxim);
//    
//#ifdef DEVELOP
//    Array *largest = Make_Array(PLAIN_KIND,UINT8_TYPE,2,maxim->dims);
//    uint8 *big     = AUINT8(largest);
//    
//    //Array_Op_Scalar(largest,SET_OP,IVAL,ValI(0));
//    
//    Array *dist    = Make_Array_With_Shape(PLAIN_KIND, UINT32_TYPE, Coord3(img.GetGridX(), img.GetGridY(), img.GetDepth()));
//    uint32 *dis    = AUINT32(dist);
//#endif
//    
//    Use_Extend_Boundary();
//    
//    Frame * fid = Make_Frame(index,Coord2(imgPara.radius,imgPara.radius),Coord2(0,0));
//    Frame * fmx = Make_Frame(maxim,Coord2(imgPara.radius,imgPara.radius),Coord2(0,0));
//    Offs_Type * off = Frame_Offsets(fid);
//    Size_Type fsz   = AForm_Size(fid);
//    
//    Histogram * wnh = Make_Histogram(UVAL,256,ValU(1),ValU(0));
//    
//    Size_Type * hcnt = wnh->counts;
//    Size_Type * dcnt = (Size_Type *) Guarded_Malloc(sizeof(Size_Type)*img.GetDepth(), Program_Name());
//    
//    Indx_Type lp = 0;
//    for (Dimn_Type y = 0; y <= img.GetHeight()-imgPara.radius; y += imgPara.radius){
//      for (Dimn_Type x = 0; x <= img.GetWidth()-imgPara.radius; x += imgPara.radius){
//        int thr;
//        
//        Indx_Type p = y*img.GetWidth() + x;
//        
//        Place_Frame(fmx,p);
//        Place_Frame(fid,p);
//        
//        Empty_Histogram(wnh);
//        Histagain_Array(wnh,fmx,0);
//        
//        { int sum = 0;
//          for (thr = 255; thr > 0; thr--)
//          { sum += hcnt[thr];
//            if (sum >= imgPara.threshold)
//              break;
//          }
//        }
//        
//        for (int d = 0; d < img.GetDepth(); d++){
//          dcnt[d] = 0;
//        }
//        
//        if (Frame_Within_Array(fid)){
//          uint8  *mx = max + p;
//          uint32 *ix = idx + p;
//#ifdef DEVELOP
//          uint8  *bg = big + p;
//#endif
//          for (int d = 0; d < fsz; d++)
//            if (mx[off[d]] >= thr){
//              dcnt[ix[off[d]]] += 1;
//#ifdef DEVELOP
//              bg[off[d]] = mx[off[d]];
//#endif
//            }
//        } else {
//          std::cout << "Frame should never be outside array" << std::endl;
//          exit (0);
//        }
//        
//        { int best = 0;
//          int plan;
//          
//          for (int d = 0; d < img.GetDepth(); d++){
//            if (dcnt[d] > best){
//              best = dcnt[d];
//              plan = d;
//            }
//          }
//          lp = Coord2IdxA(level, Coord2(y/imgPara.radius, x/imgPara.radius));
//          lev[lp] = plan;
//          
//#ifdef DEVELOP
//          Indx_Type pDist;
//          
//          for (int j=0; j<img.GetDepth(); j++) {
//            pDist = Coord2IdxA(dist, Coord3(x/imgPara.radius, y/imgPara.radius, j));
//            dis[pDist] = (uint8)dcnt[j];
//          }
//#endif
//        }
//        
//      }
//    }
//    
//    free(dcnt);
//    Kill_Histogram(wnh);
//    Kill_Frame(fid);
//    Kill_Frame(fmx);
//    
//#ifdef DEVELOP
//    ShowArray(img, imgPara, "largest.tif",largest);
//    Free_Array(largest);
//    
//    Write_Image("dist.tif", dist, DONT_PRESS);
//    Free_Array(dist);
//#endif
//    
//    if (imgPara.verbose){
//      std::cout << "Z-planes selected" << std::endl;
//    }
//    
//    return level;
//  }
//  
//  /***********************
//   * Select the z-plane with the most occurrences of bright pixel
//   *    Uses a window size twice the size of the input radius
//   **********************/
//  
//  
//  Array * PlaneSelection_DoubleRadius (Image & img, ParaSet & imgPara, Array *maxim, Array *index, int rad){
//    
//    int tmpGridX = ceil ((double)img.GetWidth()/(double)(rad/2));
//    int tmpGridY  = ceil ((double)img.GetHeight()/(double)(rad/2));
//    
//    Array  *level = Make_Array_With_Shape(PLAIN_KIND,UINT32_TYPE,Coord2(tmpGridY,tmpGridX));
//    uint32 *lev   = AUINT32(level);
//    
//    uint32 *idx   = AUINT32(index);
//    uint8  *max   = AUINT8(maxim);
//    
//#ifdef DEVELOP
//    Array *largest = Make_Array(PLAIN_KIND,UINT8_TYPE,2,maxim->dims);
//    uint8 *big     = AUINT8(largest);
//    
//    //Array_Op_Scalar(largest,SET_OP,IVAL,ValI(0));
//    
//    Array *dist    = Make_Array_With_Shape(PLAIN_KIND, UINT32_TYPE, Coord3(ceil ((double)img.GetHeight()/(double)(rad)), ceil ((double)img.GetWidth()/(double)(rad)), img.GetDepth()));
//    uint32 *dis    = AUINT32(dist);
//#endif
//    
//    Use_Extend_Boundary();
//    
//    Frame * fid = Make_Frame(index,Coord2(rad,rad),Coord2(0,0));
//    Frame * fmx = Make_Frame(maxim,Coord2(rad,rad),Coord2(0,0));
//    Offs_Type * off = Frame_Offsets(fid);
//    Size_Type fsz = AForm_Size(fid);
//    
//    Histogram * wnh = Make_Histogram(UVAL,256,ValU(1),ValU(0));
//    
//    Size_Type * hcnt = wnh->counts;
//    Size_Type * dcnt = (Size_Type *) Guarded_Malloc(sizeof(Size_Type)*img.GetDepth(),Program_Name());
//    
//    
//    Indx_Type lp = 0;
//    for (Dimn_Type y = 0; y <= img.GetHeight()-rad; y += rad){
//      for (Dimn_Type x = 0; x <= img.GetWidth()-rad; x += rad){
//        int thr;
//        
//        Indx_Type p = y*img.GetWidth() + x;
//        
//        Place_Frame(fmx,p);
//        Place_Frame(fid,p);
//        
//        Empty_Histogram(wnh);
//        Histagain_Array(wnh,fmx,0);
//        
//        { int sum = 0;
//          for (thr = 255; thr > 0; thr--){
//            sum += hcnt[thr];
//            if (sum >= imgPara.threshold){
//              break;
//            }
//          }
//        }
//        
//        for (int d = 0; d < img.GetDepth(); d++){
//          dcnt[d] = 0;
//        }
//        
//        if (Frame_Within_Array(fid)){
//          uint8  *mx = max + p;
//          uint32 *ix = idx + p;
//#ifdef DEVELOP
//          uint8  *bg = big + p;
//#endif
//          for (int d = 0; d < fsz; d++)
//            if (mx[off[d]] >= thr){
//              dcnt[ix[off[d]]] += 1;
//#ifdef DEVELOP
//              bg[off[d]] = mx[off[d]];
//#endif
//            }
//        } else {
//          std::cout << "Frame should never be outside array" << std::endl;
//          exit (0);
//        }
//        
//        { int plan;
//          int best = 0;
//          
//          for (int d = 0; d < img.GetDepth(); d++){
//            if (dcnt[d] > best){
//              best = dcnt[d];
//              plan = d;
//            }
//          }
//          lp = Coord2IdxA(level, Coord2((2*y)/rad, (2*x)/rad));
//          lev[lp] = plan;
//          
//          lp = Coord2IdxA(level, Coord2((2*y+rad)/rad, (2*x)/rad));
//          lev[lp] = plan;
//          
//          lp = Coord2IdxA(level, Coord2((2*y)/rad, (2*x+rad)/rad));
//          lev[lp] = plan;
//          
//          lp = Coord2IdxA(level, Coord2((2*y+rad)/rad, (2*x+rad)/rad));
//          lev[lp] = plan;
//          
//          
//#ifdef DEVELOP
//          Indx_Type pDist;
//          
//          for (int j=0; j<img.GetDepth(); j++) {
//            pDist = Coord2IdxA(dist, Coord3(y/rad, x/rad, j));
//            dis[pDist] = (uint8)dcnt[j];
//          }
//#endif
//        }
//        
//      }
//    }
//    
//    free(dcnt);
//    Kill_Histogram(wnh);
//    Kill_Frame(fid);
//    Kill_Frame(fmx);
//    
//#ifdef DEVELOP
//    ShowArray(img, imgPara, "largest.tif",largest);
//    Free_Array(largest);
//    
//    Write_Image("dist.tif", dist, DONT_PRESS);
//    Free_Array(dist);
//    
//    ShowArray(img, imgPara, "level_inital.tif", level);
//#endif
//    
//    if (imgPara.verbose){
//      std::cout << "Z-planes selected" << std::endl;
//    }
//    
//    return level;
//  }
//  
//  /***********************
//   * Select the z-plane with the highest entropy
//   *    Uses a subset of the stack defined by a reference height map and a distance
//   **********************/
//  
//  
//  Array * SelectPlanes_Variance(Image & img, ParaSet & imgPara, Array * substack){
//    
//    Array  *level = Make_Array_With_Shape(PLAIN_KIND,UINT32_TYPE,Coord2(img.GetGridY(),img.GetGridX()));
//    uint32 *lev   = AUINT32(level);
//    
//    Use_Reflective_Boundary();
//    
//    Frame * f = Make_Frame(substack,Coord3(1,imgPara.radius,imgPara.radius),Coord3(0, 0, 0));
//    Histogram * h = Make_Histogram(UVAL,256,ValU(1),ValU(0));
//    
//    Place_Frame(f, 0);
//    
//    double maxVar;
//    int best;
//    Indx_Type p;
//    
//    for (Dimn_Type y = 0; y <= img.GetHeight()-imgPara.radius; y += imgPara.radius){
//      for (Dimn_Type x = 0; x <= img.GetWidth()-imgPara.radius; x += imgPara.radius){
//        best = 0;
//        maxVar = 0;
//        
//        for (Dimn_Type z=0; z<img.GetDepth(); z++) {
//          
//          p = Coord2IdxA(substack, Coord3(z, y, x));
//          Place_Frame(f, p);
//          
//          Empty_Histogram(h);
//          Histagain_Array(h,f,0);
//          
//          if (maxVar < Histogram_Variance(h)) {
//            maxVar = Histogram_Variance(h);
//            best = z;
//          }
//          
//        }
//        p = Coord2IdxA(level, Coord2(y/imgPara.radius, x/imgPara.radius));
//        lev[p] = best;
//        
//      }
//    }
//    
//    
//#ifdef DEVELOP
//    ShowArray(img, imgPara, "level_Variance.tif", level);
//#endif
//    if (imgPara.verbose){
//      std::cout << "Optimal planes selected (showing highest variance)" << std::endl;
//    }
//    
//    Free_Frame(f);
//    Free_Histogram(h);
//    return level;
//  }
//  
//  /***********************
//   *
//   * Select the z-plane with the highest entropy
//   *    Uses a subset of the stack defined by a reference height map and a distance
//   * 16-bit images
//   **********************/
//  
//  
//  Array * SelectPlanes_Variance_16bit(Image & img, ParaSet & imgPara, Array * substack){
//    
//    Array  *level = Make_Array_With_Shape(PLAIN_KIND,UINT32_TYPE,Coord2(img.GetGridY(),img.GetGridX()));
//    uint32 *lev   = AUINT32(level);
//    
//    Use_Reflective_Boundary();
//    
//    Frame * f = Make_Frame(substack,Coord3(1,imgPara.radius,imgPara.radius),Coord3(0, 0, 0));
//    
//    /*****
//     * Fast version
//     *****/
//    Histogram * h = Make_Histogram(UVAL,256,ValU(1),ValU(0));
//    
//    /*****
//     * More accurate but slow version
//     *****/
//    //Histogram * h = Make_Histogram(UVAL,2048,ValU(32),ValU(0));
//    
//    
//    Place_Frame(f, 0);
//    
//    double maxVar;
//    int best;
//    Indx_Type p;
//    
//    for (Dimn_Type y = 0; y <= img.GetHeight()-imgPara.radius; y += imgPara.radius){
//      for (Dimn_Type x = 0; x <= img.GetWidth()-imgPara.radius; x += imgPara.radius){
//        best = 0;
//        maxVar = 0;
//        
//        for (Dimn_Type z=0; z<img.GetDepth(); z++) {
//          p = Coord2IdxA(substack, Coord3(z, y, x));
//          Place_Frame(f, p);
//          
//          Empty_Histogram(h);
//          Histagain_Array(h,f,0);
//          
//          if (maxVar < Histogram_Variance(h)) {
//            maxVar = Histogram_Variance(h);
//            best = z;
//          }
//          
//        }
//        p = Coord2IdxA(level, Coord2(y/imgPara.radius, x/imgPara.radius));
//        lev[p] = best;
//        
//      }
//    }
//    
//    
//#ifdef DEVELOP
//    ShowArray(img, imgPara, "level_Variance.tif", level);
//#endif
//    if (imgPara.verbose){
//      std::cout << "Optimal planes selected (showing highest variance)" << std::endl;
//    }
//    
//    Free_Frame(f);
//    Free_Histogram(h);
//    return level;
//  }
//  
//  
//  /***********************
//   * Select the z-plane with the highest entropy
//   *    Uses a subset of the stack defined by a reference height map and a distance
//   **********************/
//  
//  
//  Array * SelectPlanes_Entropy(Image & img, ParaSet & imgPara, Array * substack){
//    
//    Array  *level = Make_Array_With_Shape(PLAIN_KIND,UINT32_TYPE,Coord2(img.GetGridY(),img.GetGridX()));
//    uint32 *lev   = AUINT32(level);
//    
//    Use_Reflective_Boundary();
//    
//    Frame * f = Make_Frame(substack,Coord3(1,imgPara.radius,imgPara.radius),Coord3(0, 0, 0));
//    Histogram * h = Make_Histogram(UVAL,256,ValU(1),ValU(0));
//    
//    Place_Frame(f, 0);
//    
//    double maxVar;
//    int best;
//    Indx_Type p;
//    
//    for (Dimn_Type y = 0; y <= img.GetHeight()-imgPara.radius; y += imgPara.radius){
//      for (Dimn_Type x = 0; x <= img.GetWidth()-imgPara.radius; x += imgPara.radius){
//        best = 0;
//        maxVar = 0;
//        
//        for (Dimn_Type z=0; z<img.GetDepth(); z++) {
//          p = Coord2IdxA(substack, Coord3(z, y, x));
//          Place_Frame(f, p);
//          
//          Empty_Histogram(h);
//          Histagain_Array(h,f,0);
//          
//          if (maxVar < Histogram_Entropy(h)) {
//            maxVar = Histogram_Entropy(h);
//            best = z;
//          }
//          
//        }
//        p = Coord2IdxA(level, Coord2(y/imgPara.radius, x/imgPara.radius));
//        lev[p] = best;
//        
//      }
//    }
//    
//    
//#ifdef DEVELOP
//    ShowArray(img, imgPara, "level_Entropy.tif", level);
//#endif
//    if (imgPara.verbose){
//      std::cout << "Optimal planes selected (showing highest entropy)" << std::endl;
//    }
//    
//    Free_Frame(f);
//    Free_Histogram(h);
//    return level;
//  }
//  
//  
//  /***********************
//   * Select the z-plane with the highest LOG score
//   **********************/
//  
//  
//  Array * SelectPlanes_LOG(Image & img, ParaSet & imgPara, Array * substack){
//    
//    Convert_Array_Inplace(substack, PLAIN_KIND, INT32_TYPE, 32, 0);
//    
//    Array * logImage = LOG_Array(substack, 1.0, 3);
//    
//    Convert_Array_Inplace(logImage, PLAIN_KIND, UINT32_TYPE, 32, 0);
//    
//    uint32 * logVal = AUINT32(logImage);
//    
//    
//    // Going through the image, pixel by pixel and determine the maximal absolut value for each z-column
//    
//    Array * indexLoG = Make_Array_With_Shape(PLAIN_KIND, UINT32_TYPE, Coord2(img.GetHeight(), img.GetWidth()));
//    uint32 * idxLoG = AUINT32(indexLoG);
//    
//    Array * maximLoG = Make_Array_With_Shape(PLAIN_KIND, UINT32_TYPE, Coord2(img.GetHeight(), img.GetWidth()));
//    uint32 * maxLoG = AUINT32(maximLoG);
//    
//    uint32 maxValue, maxIndex;
//    Indx_Type q;
//    
//    for (Indx_Type p = 0; p < img.GetArea(); p++){
//      maxValue = abs(logVal[p]);
//      maxIndex = 0;
//      q = p + img.GetArea();
//      
//      for (Dimn_Type d = 1; d < img.GetDepth(); d++){
//        if (abs(logVal[q]) >= maxValue) {
//          maxValue = abs(logVal[q]);
//          maxIndex = d;
//        }
//        q += img.GetArea();
//      }
//      maxLoG[p] = maxValue;
//      idxLoG[p] = maxIndex;
//    }
//    
//    Image imgLoG = Image(logImage);
//    imgLoG.SetGrid(imgPara.radius);
//    
//    Array * level = PlaneSelection(imgLoG, imgPara, maximLoG, indexLoG);
//    
//    
//#ifdef DEVELOP
//    Write_Image("logImage.tif", logImage, DONT_PRESS);
//    
//    ShowArray(imgLoG, imgPara, "maxim_loG.tif", maximLoG);
//    ShowArray(imgLoG, imgPara, "index_loG.tif", indexLoG);
//    ShowArray(img, imgPara, "level_LoG.tif", level);
//#endif
//    
//    if (imgPara.verbose) {
//      std::cout << "Optimal planes selected (using LoG)" << std::endl;
//    }
//    
//    Free_Array(maximLoG);
//    Free_Array(indexLoG);
//    Free_Array(logImage);
//    
//    return level;
//  }
//  
//  
//  /***********************
//   * Select the z-plane with the highest LOG score -> basal cell layer
//   **********************/
//  
//  
//  Array * SelectPlanes_LOG_basal(Image & img, ParaSet & imgPara, Array * substack){
//    
//    Convert_Array_Inplace(substack, PLAIN_KIND, INT32_TYPE, 32, 0);
//    
//    Array * logImage = LOG_Array(substack, 1.0, 3);
//    
//    Convert_Array_Inplace(logImage, PLAIN_KIND, UINT32_TYPE, 32, 0);
//    
//    uint32 * logVal = AUINT32(logImage);
//    
//    Write_Image("logImage.tif", logImage, DONT_PRESS);
//    
//    
//    // Going through the image, pixel by pixel and determine the maximal absolut value for each z-column
//    
//    Array * indexLoG = Make_Array_With_Shape(PLAIN_KIND, UINT32_TYPE, Coord2(img.GetHeight(), img.GetWidth()));
//    uint32 * idxLoG = AUINT32(indexLoG);
//    
//    Array * maximLoG = Make_Array_With_Shape(PLAIN_KIND, UINT32_TYPE, Coord2(img.GetHeight(), img.GetWidth()));
//    uint32 * maxLoG = AUINT32(maximLoG);
//    
//    uint32 maxValue, maxIndex;
//    Indx_Type q;
//    
//    for (Indx_Type p = 0; p < img.GetArea(); p++){
//      maxValue = abs(logVal[p]);
//      maxIndex = 0;
//      q = p + img.GetArea();
//      
//      for (Dimn_Type d = 1; d < img.GetDepth(); d++){
//        if (abs(logVal[q]) >= maxValue) {
//          maxValue = abs(logVal[q]);
//          maxIndex = d;
//        }
//        
//        q += img.GetArea();
//      }
//      maxLoG[p] = maxValue;
//      idxLoG[p] = maxIndex;
//    }
//    Image imgLoG = Image(logImage);
//    
//    
//#ifdef DEVELOP
//    ShowArray(imgLoG, imgPara, "maxim_loG.tif", maximLoG);
//    ShowArray(imgLoG, imgPara, "index_loG.tif", indexLoG);
//#endif
//    
//    
//    Free_Array(maximLoG);
//    
//    Free_Array(logImage);
//    
//    
//    // Determining the level with the highest frequency of pixels
//    
//    Array * level = Make_Array_With_Shape(PLAIN_KIND, UINT32_TYPE, Coord2(img.GetGridY(), img.GetGridX()));
//    uint32 * lev = AUINT32(level);
//    
//    Use_Reflective_Boundary();
//    
//    Frame * f = Make_Frame(indexLoG, Coord2(imgPara.radius,imgPara.radius),Coord2(0, 0));
//    Histogram * h = Make_Histogram(UVAL, img.GetDepth(), ValU(1), ValU(0));
//    Place_Frame(f, 0);
//    
//    for (Dimn_Type y = 0; y <= img.GetHeight()-imgPara.radius; y += imgPara.radius){
//      for (Dimn_Type x = 0; x <= img.GetWidth()-imgPara.radius; x += imgPara.radius){
//        
//        Place_Frame(f, Coord2IdxA(indexLoG, Coord2(y, x)));
//        Empty_Histogram(h);
//        Histagain_Array(h, f, 0);
//        
//        Size_Type * hValues = h->counts;
//        
//        int32 maxVal = 0;
//        int32 maxPlane = 0;
//        
//        for (int i = 0; i < img.GetDepth(); i++) {
//          if (hValues[i] > maxVal) {
//            maxVal = hValues[i];
//            maxPlane = (uint32) i;
//          }
//        }
//        
//        Indx_Type p = Coord2IdxA(level, Coord2(y/imgPara.radius, x/imgPara.radius));
//        lev[p] = maxPlane;
//      }
//    }
//    
//    
//#ifdef DEVELOP
//    ShowArray(img, imgPara, "level_loG.tif", level);
//#endif
//    
//    if (imgPara.verbose) {
//      std::cout << "Optimal planes selected (using LoG)" << std::endl;
//    }
//    
//    
//    Free_Frame(f);
//    Free_Histogram(h);
//    Free_Array(indexLoG);
//    
//    return level;
//  }
//  
//}