//
////
////  Interpolation.cpp
////  proj_temp_cpp
////
////  Created by Corinna Blasse on 27.07.12.
////  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
////
//
//#include <iostream>
//
//#include "Interpolation.h"
//
//using namespace ImageClass;
//using namespace Parameter;
//
//namespace Interpolation{
//    
//    /***********************
//     * Interpolate values of median to 'corners' of its pixels in an array corners
//     *    each of whose dimensions are one more than those of median
//     *
//     **********************/
//    
//    Array *InterpolateCorners(Image & img, ParaSet & imgPara, Array *median){
//        
//        Indx_Type p, q;
//        
//        int plusX = img.GetGridX() + 1;
//        int plusY = img.GetGridY() + 1;
//        
//        Array  *corners = Make_Array_With_Shape(PLAIN_KIND,UINT32_TYPE,Coord2(plusY, plusX));
//        uint32 *map     = AUINT32(corners);
//        uint32 *med     = AUINT32(median);
//        
//        for (Dimn_Type y = 1; y < img.GetGridY(); y++){
//            p = y*plusX + 1;
//            q = y*img.GetGridX() + 1;
//            for (Dimn_Type x = 1; x < img.GetGridX(); x++, q++){
//                map[p++] = med[q] + med[q-1] + med[q-img.GetGridX()] + med[q-plusX];
//            }
//        }
//        
//        p = img.GetGridY()*plusX+1;
//        q = (img.GetGridY()-1)*img.GetGridX()+1;
//        
//        for (Dimn_Type x = 1; x < img.GetGridX(); x++, q++){
//            map[x] = (med[x] + med[x-1]) << 1;
//            map[p++] = (med[q] + med[q-1]) << 1;
//        }
//        
//        p = plusX;
//        q = img.GetGridX();
//        
//        for (Dimn_Type y = 1; y < img.GetGridY(); y++){
//            map[p] = (med[q] + med[q-img.GetGridX()]) << 1;
//            map[p+img.GetGridX()] = (med[q+img.GetGridX()-1] + med[q-1]) << 1;
//            p += plusX;
//            q += img.GetGridX();
//        }
//        
//        map[0] = med[0] << 2;
//        map[img.GetGridY()*plusX] = med[(img.GetGridY()-1)*img.GetGridX()] << 2;
//        map[img.GetGridX()] = med[img.GetGridX()-1] << 2;
//        map[img.GetGridY()*plusX+img.GetGridX()] = med[img.GetGridY()*img.GetGridX()-1] << 2;
//        
//#ifdef DEVELOP
//        ShowArray(img, imgPara, "corners.tif",corners);
//#endif
//        
//        if (imgPara.verbose){
//            std::cout << "Plane corners ready" << std::endl;
//        }
//        
//        return (corners);
//    }
//    
//    
//    /***********************
//     * Use the z-plane map corners to project image onto result, where one linearly
//     *    interprets the points between each window center.
//     * 8-bit images
//     **********************/
//    
//    Array *InterpolatePlanes_8bit(Image & img, ParaSet & imgPara, Array *maxim, Array *corners){
//        
//        Array  *result = Make_Array(PLAIN_KIND,UINT8_TYPE,2, img.GetImage()->dims);
//        uint8  *res    = AUINT8(result);
//        uint32 *map    = AUINT32(corners);
//        uint8  *val    = AUINT8(img.GetImage());
//        uint8  *wal    = val + img.GetArea();
//        
//        int plusX = img.GetGridX() + 1;
//        
//#ifdef DEVELOP
//        Array *marked  = Convert_Array(img.GetImage(),RGB_KIND,UINT8_TYPE,8,0);
//        Array *residue = Make_Array(PLAIN_KIND,UINT8_TYPE,2,img.GetImage()->dims);
//        
//        uint8 *grn = AUINT8(marked) + img.GetDepth()*img.GetArea();
//        uint8 *blu = grn + img.GetDepth()*img.GetArea();
//        uint8 *rem = AUINT8(residue);
//        uint8 *max = AUINT8(maxim);
//#endif
//        
//        int n = imgPara.radius-1;
//        
//        for (Dimn_Type y = 0; y < img.GetGridY(); y++){
//            
//            uint32 *dep = map + y*plusX;
//            
//            for (Dimn_Type x = 0; x < img.GetGridX(); x++){
//                double  ul = dep[0] / 4.;
//                double  ur = dep[1] / 4.;
//                double  ll = dep[plusX] / 4.;
//                double  lr = dep[plusX+1] / 4.;
//                
//                for (int r = 0; r <= n; r++){
//                    double ru = ul*(n-r)/n + ur*r/n;
//                    double rl = ll*(n-r)/n + lr*r/n;
//                    
//                    for (int s = 0; s <= n; s++){
//                        double    wt = ru*(n-s)/n + rl*s/n;
//                        int       wi = (int) wt;
//                        double    wf = wt-wi;
//                        Indx_Type p = (y*imgPara.radius+s)*img.GetWidth() + (x*imgPara.radius+r);
//                        Indx_Type q = (Indx_Type) wi*img.GetArea() + p;
//                        
//                        if (wi+1 < img.GetDepth())
//                            res[p] = val[q] * (1.-wf) + wal[q] * wf;
//                        else
//                            res[p] = val[q];
//#ifdef DEVELOP
//                        rem[p] = max[p] - res[p];
//                        if (wf < 1.)
//                        { grn[q] *= (1.-wf)*(1.-wf);
//                            blu[q] *= wf*wf;
//                        }
//                        if (wf > 0.)
//                        { grn[q+img.GetArea()] *= wf*wf;
//                            blu[q+img.GetArea()] *= (1.-wf)*(1.-wf);
//                        }
//#endif
//                    }
//                }
//                dep += 1;
//            }
//        }
//        
//#ifdef DEVELOP
//        ShowArray(img, imgPara, "marked.tif",marked);
//        ShowArray(img, imgPara, "residue.tif",residue);
//        Free_Array(residue);
//        Free_Array(marked);
//#endif
//        
//        if (imgPara.verbose){
//            std::cout << "Result interpolated" << std::endl;
//        }
//        
//        return (result);
//    }
//    
//    /***********************
//     * Use the z-plane map corners to project image onto result, where one linearly
//     *    interprets the points between each window center.
//     * 16-bit images
//     **********************/
//    
//    Array *InterpolatePlanes_16bit(Image & img, ParaSet & imgPara, Array *maxim, Array *corners){
//        
//        Array  *result = Make_Array(PLAIN_KIND,UINT16_TYPE,2, img.GetImage()->dims);
//        uint16  *res    = AUINT16(result);
//        uint32 *map    = AUINT32(corners);
//        uint16  *val    = AUINT16(img.GetImage());
//        uint16  *wal    = val + img.GetArea();
//        
//        int plusX = img.GetGridX() + 1;
//        
//#ifdef DEVELOP
//        //Array *marked  = Convert_Array(img.GetImage(),RGB_KIND,UINT8_TYPE,8,0);
//        Array *residue = Make_Array(PLAIN_KIND,UINT16_TYPE,2,img.GetImage()->dims);
//        
//        //uint16 *grn = AUINT16(marked) + img.GetDepth()*img.GetArea();
//        //uint16 *blu = grn + img.GetDepth()*img.GetArea();
//        uint16 *rem = AUINT16(residue);
//        uint16 *max = AUINT16(maxim);
//#endif
//        
//        int n = imgPara.radius-1;
//        
//        for (Dimn_Type y = 0; y < img.GetGridY(); y++){
//            
//            uint32 *dep = map + y*plusX;
//            
//            for (Dimn_Type x = 0; x < img.GetGridX(); x++){
//                double  ul = dep[0] / 4.;
//                double  ur = dep[1] / 4.;
//                double  ll = dep[plusX] / 4.;
//                double  lr = dep[plusX+1] / 4.;
//                
//                for (int r = 0; r <= n; r++){
//                    double ru = ul*(n-r)/n + ur*r/n;
//                    double rl = ll*(n-r)/n + lr*r/n;
//                    
//                    for (int s = 0; s <= n; s++){
//                        double    wt = ru*(n-s)/n + rl*s/n;
//                        int       wi = (int) wt;
//                        double    wf = wt-wi;
//                        Indx_Type p = (y*imgPara.radius+s)*img.GetWidth() + (x*imgPara.radius+r);
//                        Indx_Type q = (Indx_Type) wi*img.GetArea() + p;
//                        
//                        if (wi+1 < img.GetDepth())
//                            res[p] = val[q] * (1.-wf) + wal[q] * wf;
//                        else
//                            res[p] = val[q];
//#ifdef DEVELOP
//                        rem[p] = max[p] - res[p];
//                        //if (wf < 1.)
//                        //{ grn[q] *= (1.-wf)*(1.-wf);
//                        //  blu[q] *= wf*wf;
//                        //}
//                        //if (wf > 0.)
//                        //{ grn[q+img.GetArea()] *= wf*wf;
//                        //  blu[q+img.GetArea()] *= (1.-wf)*(1.-wf);
//                        //}
//#endif
//                    }
//                }
//                dep += 1;
//            }
//        }
//        
//#ifdef DEVELOP
//        //ShowArray(img, imgPara, "marked.tif",marked);
//        ShowArray(img, imgPara, "residue.tif",residue);
//        Free_Array(residue);
//        //Free_Array(marked);
//#endif
//        
//        if (imgPara.verbose){
//            std::cout << "Result interpolated" << std::endl;
//        }
//        
//        return (result);
//    }
//    
//    /***********************
//     * Use the z-plane map corners to project image onto result, where one linearly
//     *    interprets the points between each window center.
//     * 8-bit images
//     **********************/
//    
//    Array *InterpolatePlanes_8bit_MaxInterpolation(Image & img, ParaSet & imgPara, Array *maxim, Array *corners){
//        
//        Array  *result = Make_Array(PLAIN_KIND,UINT8_TYPE,2, img.GetImage()->dims);
//        uint8  *res    = AUINT8(result);
//        uint32 *map    = AUINT32(corners);
//        uint8  *val    = AUINT8(img.GetImage());
//        uint8  *wal    = val + img.GetArea();
//        
//        int plusX = img.GetGridX() + 1;
//        
//        Array *maxInt  = Convert_Array(maxim, RGB_KIND, UINT8_TYPE, maxim->scale, 0);
//        
//#ifdef DEVELOP
//        Array *marked  = Convert_Array(img.GetImage(),RGB_KIND,UINT8_TYPE,8,0);
//        Array *residue = Make_Array(PLAIN_KIND,UINT8_TYPE,2,img.GetImage()->dims);
//        
//        uint8 *grn = AUINT8(marked) + img.GetDepth()*img.GetArea();
//        uint8 *blu = grn + img.GetDepth()*img.GetArea();
//        uint8 *rem = AUINT8(residue);
//        uint8 *max = AUINT8(maxim);
//#endif
//        
//        int n = imgPara.radius-1;
//        
//        for (Dimn_Type y = 0; y < img.GetGridY(); y++){
//            
//            uint32 *dep = map + y*plusX;
//            
//            for (Dimn_Type x = 0; x < img.GetGridX(); x++){
//                double  ul = dep[0] / 4.;
//                double  ur = dep[1] / 4.;
//                double  ll = dep[plusX] / 4.;
//                double  lr = dep[plusX+1] / 4.;
//                
//                for (int r = 0; r <= n; r++){
//                    double ru = ul*(n-r)/n + ur*r/n;
//                    double rl = ll*(n-r)/n + lr*r/n;
//                    
//                    for (int s = 0; s <= n; s++){
//                        double    wt = ru*(n-s)/n + rl*s/n;
//                        int       wi = (int) wt;
//                        double    wf = wt-wi;
//                        Indx_Type p = (y*imgPara.radius+s)*img.GetWidth() + (x*imgPara.radius+r);
//                        
//                        Indx_Type q = (Indx_Type) wi*img.GetArea() + p;
//                        
//                        
//                        int lowerLayer = fmax(0, wi-2);
//                        int upperLayer = fmin(img.GetDepth()-1, wi+2);
//                        
//                        for (int currLayer = lowerLayer; currLayer <= upperLayer; currLayer++) {
//                            Indx_Type currQ = (Indx_Type) currLayer*img.GetArea() + p;
//                            if (val[currQ] > val[q]) {
//                                //std::cout << (int)val[currQ]  << " > " << (int) val[q] << std::endl;
//                                
//                                q = currQ;
//                                wi = currLayer;
//                                
//                                if (val[q] > 20) {
//                                    Draw_Point(maxInt, &PURPLE, Coord2(y*imgPara.radius+s, x*imgPara.radius+r ));
//                                }
//                                
//                            }
//                        }
//                        
//                        //if (wi+1 < img.GetDepth())
//                        if (wi+1 < img.GetDepth() && wf==wt-wi)
//                            res[p] = val[q] * (1.-wf) + wal[q] * wf;
//                        else
//                            res[p] = val[q];
//                        
//                        
//#ifdef DEVELOP
//                        rem[p] = max[p] - res[p];
//                        if (wf < 1.)
//                        { grn[q] *= (1.-wf)*(1.-wf);
//                            blu[q] *= wf*wf;
//                        }
//                        if (wf > 0.)
//                        { grn[q+img.GetArea()] *= wf*wf;
//                            blu[q+img.GetArea()] *= (1.-wf)*(1.-wf);
//                        }
//#endif
//                    }
//                }
//                dep += 1;
//            }
//        }
//        
//        //ShowArray(img, imgPara, "maxIntChange.tif", maxInt);
//        Free_Array(maxInt);
//        
//#ifdef DEVELOP
//        ShowArray(img, imgPara, "marked.tif",marked);
//        ShowArray(img, imgPara, "residue.tif",residue);
//        Free_Array(residue);
//        Free_Array(marked);
//#endif
//        
//        if (imgPara.verbose){
//            std::cout << "Result interpolated" << std::endl;
//        }
//        
//        return (result);
//    }
//    
//    
//    /***********************
//     * Use the z-plane map corners to project image onto result, where one linearly
//     *    interprets the points between each window center.
//     * 16-bit images
//     **********************/
//    
//    Array *InterpolatePlanes_16bit_MaxInterpolation(Image & img, ParaSet & imgPara, Array *maxim, Array *corners){
//        
//        Array  *result = Make_Array(PLAIN_KIND,UINT16_TYPE,2, img.GetImage()->dims);
//        uint16  *res    = AUINT16(result);
//        uint32 *map    = AUINT32(corners);
//        uint16  *val    = AUINT16(img.GetImage());
//        uint16  *wal    = val + img.GetArea();
//        
//        int plusX = img.GetGridX() + 1;
//        
//#ifdef DEVELOP
//        //Array *marked  = Convert_Array(img.GetImage(),RGB_KIND,UINT8_TYPE,8,0);
//        Array *residue = Make_Array(PLAIN_KIND,UINT16_TYPE,2,img.GetImage()->dims);
//        
//        //uint8 *grn = AUINT8(marked) + img.GetDepth()*img.GetArea();
//        //uint8 *blu = grn + img.GetDepth()*img.GetArea();
//        uint16 *rem = AUINT16(residue);
//        uint16 *max = AUINT16(maxim);
//#endif
//        
//        int n = imgPara.radius-1;
//        
//        for (Dimn_Type y = 0; y < img.GetGridY(); y++){
//            
//            uint32 *dep = map + y*plusX;
//            
//            for (Dimn_Type x = 0; x < img.GetGridX(); x++){
//                double  ul = dep[0] / 4.;
//                double  ur = dep[1] / 4.;
//                double  ll = dep[plusX] / 4.;
//                double  lr = dep[plusX+1] / 4.;
//                
//                for (int r = 0; r <= n; r++){
//                    double ru = ul*(n-r)/n + ur*r/n;
//                    double rl = ll*(n-r)/n + lr*r/n;
//                    
//                    for (int s = 0; s <= n; s++){
//                        double    wt = ru*(n-s)/n + rl*s/n;
//                        int       wi = (int) wt;
//                        double    wf = wt-wi;
//                        Indx_Type p = (y*imgPara.radius+s)*img.GetWidth() + (x*imgPara.radius+r);
//                        
//                        Indx_Type q = (Indx_Type) wi*img.GetArea() + p;
//                        
//                        
//                        int lowerLayer = fmax(0, wi-2);
//                        int upperLayer = fmin(img.GetDepth()-1, wi+2);
//                        
//                        for (int currLayer = lowerLayer; currLayer <= upperLayer; currLayer++) {
//                            Indx_Type currQ = (Indx_Type) currLayer*img.GetArea() + p;
//                            if (val[currQ] > val[q]) {                
//                                q = currQ;
//                                wi = currLayer;
//                                
//                            }
//                        }
//                        
//                        //if (wi+1 < img.GetDepth())
//                        if (wi+1 < img.GetDepth() && wf==wt-wi)
//                            res[p] = val[q] * (1.-wf) + wal[q] * wf;
//                        else
//                            res[p] = val[q];
//                        
//                        
//#ifdef DEVELOP
//                        rem[p] = max[p] - res[p];
//                        /*
//                         if (wf < 1.)
//                         { grn[q] *= (1.-wf)*(1.-wf);
//                         blu[q] *= wf*wf;
//                         }
//                         if (wf > 0.)
//                         { grn[q+img.GetArea()] *= wf*wf;
//                         blu[q+img.GetArea()] *= (1.-wf)*(1.-wf);
//                         }*/
//#endif
//                    }
//                }
//                dep += 1;
//            }
//        }
//        
//        
//#ifdef DEVELOP
//        //ShowArray(img, imgPara, "marked.tif",marked);
//        ShowArray(img, imgPara, "residue.tif",residue);
//        Free_Array(residue);
//        //Free_Array(marked);
//#endif
//        
//        if (imgPara.verbose){
//            std::cout << "Result interpolated" << std::endl;
//        }
//        
//        return (result);
//    }
//    
//    
//}