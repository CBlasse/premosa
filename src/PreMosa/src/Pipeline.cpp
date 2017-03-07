//
//  Pipeline.cpp
//  CMakeProject
//
//  Created by blasse on 7/26/13.
//
//

#include "Pipeline.h"

#include <algorithm>
#include <math.h>
#include <vector>
#include <dirent.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <sys/stat.h>

extern "C" {
#include "array.h"
#include "filters.h"
#include "histogram.h"
#include "image.h"
}

// Includes from HelperFunctions library
#include "MylibValues.h"
// Includes from ProjectionLibrary
#include "Parameter.h"
#include "ProjectionAlgorithm.h"
#include "Filter.h"
#include "RasterizedImage.h"

#include "FileHandling.h"

namespace PreprocessingPipeline {
  
  namespace {
    ////////////////////////////////////////////////////////////////////////////////
    /*
     Function to check whether the file exists
     
     @param filename:   filename
     */
    inline bool FileExists (const std::string& filename) {
      
      if (std::ifstream(filename.c_str()))
        {
        return true;
        }
      return false;
    }

    
    ////////////////////////////////////////////////////////////////////////////////
    /*
     Flat field correction
     -> corrected image = (image / flat field) * rescaling factor
     
     @param image:            input image (8-bit or 16-bit)
     @param flatField:        flat field image  (32-bit)
     @param rescalingFactor:  factor to shift back the corrected image to its original intensity range
     */
    
    Array * FlatFieldCorrection (Array * image, Array * flatField, double rescalingFactor) {
      
      
      if (image->type == UINT8_TYPE) {
        
        uint8 * imgVals = AUINT8(image);
        float32 * noiseVals = AFLOAT32(flatField);
        Array * img8 = Make_Array(PLAIN_KIND, UINT8_TYPE, 2, image->dims);
        uint8 * img8Vals = AUINT8(img8);
        
        for (int p=0; p<image->dims[0]*image->dims[1]; ++p) {
          double val = rescalingFactor * ((double)imgVals[p] / (double)noiseVals[p]);     // normalization
          if (val > 255) {                                                                // clipping of oversaturated values
            img8Vals[p] = 255;
          } else {
            img8Vals[p] = (uint8) val;
          }
        }
        
        return img8;
        
      } else if (image->type == UINT16_TYPE) {
        
        uint16 * imgVals = AUINT16(image);
        float32 * noiseVals = AFLOAT32(flatField);
        Array * img16 = Make_Array(PLAIN_KIND, UINT16_TYPE, 2, image->dims);
        uint16 * img16Vals = AUINT16(img16);
        
        for (int p=0; p<image->dims[0]*image->dims[1]; ++p) {
          double val = rescalingFactor * ((double)imgVals[p] / (double)noiseVals[p]);   // normalization
          if (val > 65535) {                                                            // clipping of oversaturated values
            img16Vals[p] = 65535;
          } else {
            img16Vals[p] = (uint16) val;
          }
        }
        
        return img16;
        
      } else {
        std::cout << "Flat field correction is not supported for that data type " << std::endl;
      }
      
      return image;
    }
    
    
    
    ////////////////////////////////////////////////////////////////////////////////
    /*
     Approximation of the flat field by filtering the image
     
     @param image:            input image (8-bit or 16-bit)
     */
    
    Array * FlatFieldApproximation (Array * image) {
      
      Array * approximatedFlatField = Make_Array_With_Shape(PLAIN_KIND, FLOAT32_TYPE, Coord2(image->dims[1], image->dims[0]));
      Range_Bundle range;
      Array_Range(&range, image);
      
      Use_Reflective_Boundary();
      
      // frames defines the local neighborhood
      int filterRadius = image->dims[0]/5;
      if (filterRadius%2 == 0) {
        filterRadius++;
      }
      Frame * f = Make_Frame(image,Coord2(filterRadius,filterRadius),Coord2((filterRadius-1)/2,(filterRadius-1)/2));
      Histogram * h = Make_Histogram(UVAL,range.maxval.uval,MylibValues::ValU(1), MylibValues::ValU(0));
      Place_Frame(f,0);
      
      for (Indx_Type p = 0; p < image->size; p++){
        Empty_Histogram(h);
        Histagain_Array(h,f,0);
        Set_Array_Value(approximatedFlatField, Idx2CoordA(approximatedFlatField, p), MylibValues::ValF(Percentile2Bin(h,.5)));    // median is given at 50th percentile
        Move_Frame_Forward(f);
      }
      
      Kill_Histogram(h);
      Kill_Frame(f);
      
      return approximatedFlatField;
    }
    
    
    
    ////////////////////////////////////////////////////////////////////////////////
    /*
     Method to add leading zeros to the string
     
     @param number:    number with should be padded
     @param maxDigits: maximal number of digits
     */
    
    std::string PadString (char paddingCharacter, int number, int maxDigits )
    {
    std::stringstream paddedString;
    paddedString << std::setw( maxDigits ) << std::setfill( paddingCharacter ) << number;
    return paddedString.str();
    }
    
    
    ////////////////////////////////////////////////////////////////////////////////
    /*
     Creation of a tile configuration using a master tile configuration
     
     @param inputParamter:    struct with parameters
     @param tileConfigPath:   file name of the new tile configuration file
     @param fileName:         name template of the input images
     */
    
    
    void CreateTileConfiguration (XmlFile inputParamter, std::string tileConfigPath, std::string fileName) {
      
      std::ifstream masterFile;
      masterFile.open (inputParamter.GetMasterTileConfig());
      
      if (masterFile.is_open()) {
        std::cout << "Master tile configuration file cannot be opened! A new configuration will be computed!" << std::endl;
        return;
      }
      
      int paddingNumber = inputParamter.GetFileNameExpression().find("}") - inputParamter.GetFileNameExpression().find("{") - 1;  // determination of the number of digits of the tile indication
      
      if (FileExists(tileConfigPath)) {
        system(&("rm "+ tileConfigPath)[0]);
      }
      
      std::ofstream tileConfig;
      tileConfig.open (tileConfigPath, std::ofstream::trunc); // writing of the new tile configuration file
      if (tileConfig.is_open()) {
        
        std::string line;
        while (std::getline(masterFile, line)) {
          
          for (int i=inputParamter.GetFirstIndex(); i<=inputParamter.GetGridSizeX()*inputParamter.GetGridSizeY(); i++) {  // iteration through the number of tile
            
            std::string expression = '{' + PadString('i', 'i', paddingNumber) +'}'; // replace the current tile number
            std::string tile = PadString('0', i, paddingNumber);
            std::string tileReplacement = "%TILE" + std::to_string(i) +"%";
            std::string newFileName = fileName;
            newFileName = newFileName.replace(inputParamter.GetFileNameExpression().find("{"), paddingNumber, tile);
            
            line = line.replace(inputParamter.GetFileNameExpression().find("%"), tileReplacement.length(), newFileName);
          }
          tileConfig << line;
        }
      }
      
      masterFile.close();
      tileConfig.close();
    }
    
    
    ////////////////////////////////////////////////////////////////////////////////
    /*
     Review if the current tile presents a 'black tile' (tile with no signal)
     -> the standart deviation of the intensity distribution decided about it
     -> if it is a black tile, then the tile name is added to the blackTiles vector
     
     @param image:      input tile
     @param fileName:   file name of the input image
     @param blackTiles: vector of black tiles (will be extended if the current tile features no signals)
     */
    
    Array * Convl_8bit (Array * img, Double_Matrix * kernel) {
      
      Use_Reflective_Boundary();
      
      float64 * kernelVals = AFLOAT64(kernel);
      
      Array * filteredImg = Make_Array_With_Shape(PLAIN_KIND, UINT8_TYPE, Coord2(img->dims[1], img->dims[0]));
      uint8 * filVals = AUINT8(filteredImg);
      
      Frame * f = Make_Frame(img, Coord2(3, 3), Coord2(1, 1));
      Place_Frame(f, 0);
      
      for (Indx_Type p =0; p<filteredImg->size; p++) {
        uint8 * data = (uint8 *) Frame_Values(f);
        
        double sum = 0;
        for (int i=0; i<9; i++) {
          sum += double (data[i]) * double (kernelVals[i]);
        }
        filVals[p] = (uint8) sum;
        Move_Frame_Forward(f);
      }
      Free_Frame(f);
      
      return filteredImg;
    }
    
    Array * Convl_16bit (Array * img, Double_Matrix * kernel) {
      
      Use_Reflective_Boundary();
      
      float64 * kernelVals = AFLOAT64(kernel);
      
      Array * filteredImg = Make_Array_With_Shape(PLAIN_KIND, UINT16_TYPE, Coord2(img->dims[1], img->dims[0]));
      uint16 * filVals = AUINT16(filteredImg);
      
      Frame * f = Make_Frame(img, Coord2(3, 3), Coord2(1, 1));
      Place_Frame(f, 0);
      
      for (Indx_Type p =0; p<filteredImg->size; p++) {
        uint16 * data = (uint16 *) Frame_Values(f);
        
        double sum = 0;
        for (int i=0; i<9; i++) {
          sum += double (data[i]) * double (kernelVals[i]);
        }
        filVals[p] = (uint16) sum;
        Move_Frame_Forward(f);
      }
      Free_Frame(f);
      
      return filteredImg;
    }

    

    
    void CheckForBlackTile (Array * image, std::string fileName, std::vector<std::string> &blackTiles) {
      
      Histogram *h;
 
      Double_Array * h1 = Filter_Power(Gaussian_Filter(1,1),2);
      
      if (image->type == UINT8_TYPE) {
        Array * gauss1 = Convl_8bit(image, h1);

        h = Histogram_Array(gauss1, 256, MylibValues::ValU(0), MylibValues::ValU(1));
        if (Histogram_Sigma(h) < 2 ) {          // Decision based on a fixed threshold (~0.6% of 256)
          std::cout << fileName << std::endl;
          blackTiles.push_back(fileName);
        }
      } else {
        Array * gauss1 = Convl_16bit(image, h1);

        h = Histogram_Array(gauss1, 65536, MylibValues::ValU(0), MylibValues::ValU(1));
        if (Histogram_Sigma(h) < 150 ) {         // Decision based on a fixed threshold (~0.02% of 65536)
          std::cout << fileName << std::endl;
          blackTiles.push_back(fileName);
        }
      }
      
      Free_Histogram(h);
    }
    
    
    ////////////////////////////////////////////////////////////////////////////////
    /*
     Function to remove all black tiles from the tile configuration (necessary to enable the contrast adjustment)
     
     @param tileConfigPath:   original tile configuration path
     @param tileConfigPath2:  new tile configuration path
     @param blackTiles:       vector with the file names of all black tiles
     */
    
    void RemoveBlackTilesFromTileConfiguration (std::string tileConfigPath,std::string tileConfigPath2, std::vector<std::string> &blackTiles) {
      
      std::ifstream masterFile;
      masterFile.open (tileConfigPath);
      
      if (FileExists(tileConfigPath2)) {
        system(&("rm "+ tileConfigPath2)[0]);
      }
      
      std::ofstream masterFile2;
      masterFile2.open (tileConfigPath2);
      
      if (!masterFile.is_open() or !masterFile2.is_open()) {
        return;
      }
      
      std::string line;

      while (std::getline(masterFile, line)) {
        
        std::string fileName = line.substr(0, line.find(";"));
        if (std::find(blackTiles.begin(), blackTiles.end(), fileName) == blackTiles.end()) {    // if the file name is not included in the black tile vector, then the line is added to the new tile configuration
          masterFile2 << line << "\n";
        }
      }
      
      masterFile.close();
      masterFile2.close();
    }
    
    
    ////////////////////////////////////////////////////////////////////////////////
    /*
     Function to check whether the filename has ".tif" as ending
     
     @param file:   file
    */
    
    bool CheckFileEnding (const struct dirent *file) {
      std::string fileEnding = ".tif";
      std::string filename = file->d_name;
      
      if (filename.find(fileEnding) != std::string::npos) {
        return true;
      }
      
      return false;
    }
    
    
  } // namespace
  
  


  
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   IMAGE PREPROCESSING PIPELINE
   
   -> Subsequent processing of the individual time points
   - Projection (8bit / 16bit)
   - Flat Field Correction (8bit / 16bit)
   - Contrast Optimization (8bit / 16bit) / Fiji
   - Stitching (8bit / 16bit) / Fiji
   
   @param inputParameter:  struct with all parameters
   */
  
  void Pipeline (XmlFile inputParamter) {
    
    if (inputParamter.GetComplete()) {    // Required parameter have to be set
      
      DIR *inputDir;
      inputDir = opendir(&(inputParamter.GetInputDir()[0]));
      OutputDirs outputDirs = MakeOutputDirs (inputParamter.GetOutputDir());  // Creation of the required output folders
      
      int paddingNumber = inputParamter.GetFolderNameExpression().find("}") - inputParamter.GetFolderNameExpression().find("{") - 1;    // Parsing the number of digits used to encode the time point
      std::string expression = '{' + PadString('i', 'i', paddingNumber)+'}';
      
      
      // Settings for the flat field correction
      //  1) no flat field correction (empty string in xml)
      //  2) flat field correction using a flat field approximated from the actual image ("Approximation" in xml)
      //  3) flat field correction using a measured flat field (file name in xml)
      bool flatFieldCorrection = true;
      bool withoutFlatField = false;
      double flatFieldMaximum;
      Array * flatField;
      
      if (inputParamter.GetFlatField() == "") {
        flatFieldCorrection = false;
      } else if (inputParamter.GetFlatField() == "Approximation") {
        withoutFlatField = true;
      } else {
        flatField = Read_Image(&(inputParamter.GetFlatField()[0]),0);
        Range_Bundle flatFieldRange;
        Array_Range(&flatFieldRange, flatField);
        flatFieldMaximum = ceil(flatFieldRange.maxval.fval);
      }
      
      
      if (inputDir) {
        
        // Initialize the parameter struct
        ProjectionMethod::ParaSet projectionParameters;
        projectionParameters.verbose = (bool) inputParamter.GetVerbose();
        projectionParameters.radius = inputParamter.GetRadius();
        projectionParameters.distance = inputParamter.GetDistance();
        projectionParameters.threshold = inputParamter.GetThreshold();
        projectionParameters.printRealHeightMap = (bool) inputParamter.GetPrintHeightMap();
        int originalRadius = inputParamter.GetRadius();
        
        if (projectionParameters.radius % 2 != 0) {
          projectionParameters.radius +=1;
          originalRadius +=1;
        }
        
        // Set default minIntensity and maxIntensity values for the contrast adjustement (it depends on the data type, therefore it could not be set up before)
        if (inputParamter.GetMinIntensity() == 0 and inputParamter.GetMaxIntensity() == 65535) {
          std::string iPadded = PadString('0', inputParamter.GetTStart(), paddingNumber);
          std::string folderName (inputParamter.GetFolderNameExpression());
          folderName.replace(folderName.find("{"), expression.length(), iPadded);
          DIR * exampleDir;
          exampleDir = opendir(&(inputParamter.GetInputDir()+"/"+folderName)[0]);
          struct dirent *exampleEntry;
          bool firstFile = true;
          std::string exampleFileName;
          
          while (firstFile and (exampleEntry = readdir(exampleDir)) != NULL){
            exampleFileName = inputParamter.GetInputDir()+"/"+folderName + "/"+ exampleEntry->d_name;
          }
          
          Array * exampleImage = Read_Image(&exampleFileName[0],0);
          if (exampleImage->type == UINT8_TYPE) {
            inputParamter.SetMinIntensity(10);
            inputParamter.SetMaxIntensity(245);
          }
          if (exampleImage->type == UINT16_TYPE) {
            inputParamter.SetMinIntensity(50);
            inputParamter.SetMaxIntensity(65000);
          }
          Free_Array(exampleImage);
        }
        
        for (int i=inputParamter.GetTStart(); i<=inputParamter.GetTEnd(); ++i) {
          
          std::string iPadded = PadString('0', inputParamter.GetTStart(), paddingNumber);
          std::string folderName (inputParamter.GetFolderNameExpression());
          folderName.replace(folderName.find("{"), expression.length(), iPadded);
          
          std::cout << "Processed time point: " << &(inputParamter.GetInputDir()+"/"+folderName)[0] << std::endl;
          DIR * timePointDir;
          timePointDir = opendir(&(inputParamter.GetInputDir()+"/"+folderName)[0]);
          
          std::vector<std::string> blackTiles;
          
          if (timePointDir != NULL) {
            
            struct dirent *entry;
            
            // Create subdirectory
            DIR * outputDirFFC_TP;
            std::string outputDirFFC_TP_Path = outputDirs.flatFieldPath+"/"+folderName;
            outputDirFFC_TP = opendir(&(outputDirFFC_TP_Path)[0]);
            if (outputDirFFC_TP == NULL) {
              mkdir(&(outputDirFFC_TP_Path)[0], 0777);
              outputDirFFC_TP  = opendir (&(outputDirFFC_TP_Path)[0]);
            }
            
            while ((entry = readdir(timePointDir)) != NULL){
              
              std::string entryString = entry->d_name;
              if (entryString == "." or entryString == "..") {
                continue;
              }
              // Read the image
              std::string path = inputParamter.GetInputDir()+"/"+folderName+"/"+entry->d_name;
              std::cout << path << std::endl;
              Array * image = Read_Image(&path[0],0);
              
              // Resetting the radius parameter
              projectionParameters.radius = originalRadius;
              
              // Stack projection (if the input mosaic is only 2D, then the images are just copied)
              Array * projection;

              if (image->ndims == 2) {
                projection = Copy_Array(image);
              } else {
                ProjectionMethod::ProjectionAlgorithm projectionAlgorithm (image);
                projectionAlgorithm.SetParameters(projectionParameters);
                projectionAlgorithm.ProjectImageStack();
                projection = projectionAlgorithm.GetProjection(0);

                if (inputParamter.GetPrintHeightMap()) {
                  projectionAlgorithm.DrawInterpolatedHeightMap(&(outputDirFFC_TP_Path+"/"+entry->d_name)[0]);
                }
              }

              
              // Review if the image features significant signals and is not a 'black tile'
              CheckForBlackTile(projection, entry->d_name, blackTiles);
              
              // Flat field correction
              if (flatFieldCorrection) {
                if (withoutFlatField) {
                  Free_Array(flatField);
                  flatField = FlatFieldApproximation(projection);
                  Range_Bundle flatFieldRange;
                  Array_Range(&flatFieldRange, flatField);
                  flatFieldMaximum = ceil(flatFieldRange.maxval.fval);
                }
                Array * correctedImg = FlatFieldCorrection(projection, flatField, flatFieldMaximum);
                Write_Image(&(outputDirFFC_TP_Path+"/"+entry->d_name)[0], correctedImg, DONT_PRESS);
                Free_Array(correctedImg);
              } else {
                Write_Image(&(outputDirFFC_TP_Path+"/"+entry->d_name)[0], projection, DONT_PRESS);
              }
              
              // Cleaning (projection is freed by the destructor of the ProjectionAlgorithm)
              // Free_Array(projection);
              Free_Array(image);
            }
            
            // Folder for the contrast adjustment
            DIR * outputDirAC_TP;
            std::string outputDirAC_TP_Path = outputDirs.contrastAdjPath+"/"+folderName;
            outputDirAC_TP = opendir(&(outputDirAC_TP_Path)[0]);
            if (outputDirAC_TP == NULL) {
              mkdir(&(outputDirAC_TP_Path)[0], 0777);
              outputDirAC_TP  = opendir (&(outputDirAC_TP_Path)[0]);
            }
            
            std::string fileName (inputParamter.GetFileNameExpression());
            fileName = fileName.replace(fileName.find("TIME"), 4, iPadded);
            
            // Mosaic stitching to determine the tile configuration
            std::string tileConfig = outputDirFFC_TP_Path + "/" + "TileConfiguration.txt";
            if (inputParamter.GetMasterTileConfig() != "") {
              CreateTileConfiguration(inputParamter, tileConfig, fileName);
            }
            
            std::string stitchProgram;
            
            // Fiji call (depends on the presence of a master tile configuration and the system(Linux/OsX) )
            if (!FileExists(tileConfig)) {
              
#ifdef __APPLE__
              stitchProgram = inputParamter.GetFiji() + " -batch" +
              " " + inputParamter.GetScriptLocation() + "/gridStitching_GetTileConfig_MacOsX.bsh" +
              " -i" + outputDirFFC_TP_Path +
              "=-x" + std::to_string(inputParamter.GetGridSizeX()) +
              "=-y" + std::to_string(inputParamter.GetGridSizeY()) +
              "=-l" +std::to_string(inputParamter.GetOverlap())+
              "=-n" + fileName +
              "=-f"+ std::to_string(inputParamter.GetFirstIndex());
#endif
              
#ifdef __linux__
              stitchProgram = "xvfb-run -a " + inputParamter.GetFiji() +
              " -Di=" + outputDirFFC_TP_Path +
              " -Dx=" + std::to_string(inputParamter.GetGridSizeX()) +
              " -Dy=" + std::to_string(inputParamter.GetGridSizeY()) +
              " -Dl=" +std::to_string(inputParamter.GetOverlap())+
              " -Dn=" + fileName +
              " -Df="+ std::to_string(inputParamter.GetFirstIndex()) +
              " -- --no-splash " + inputParamter.GetScriptLocation() + "/gridStitching_GetTileConfig.bsh";
#endif
              
            } else {
#ifdef __APPLE__
              stitchProgram = inputParamter.GetFiji() + " -batch" +
              " " + inputParamter.GetScriptLocation() + "/gridStitching_UseTileConfig_MacOsX.bsh"+
              " -in" + outputDirFFC_TP_Path +
              "=-out" + outputDirFFC_TP_Path+"/";
              
#endif
              
#ifdef __linux__
              stitchProgram = "xvfb-run -a " + inputParamter.GetFiji() +
              " -Din=" + outputDirFFC_TP_Path +
              " -Dout=" + outputDirFFC_TP_Path+ "/" +
              " -- --no-splash " + inputParamter.GetScriptLocation() + "/gridStitching_UseTileConfig.bsh";
#endif
            }
  
            system(&stitchProgram[0]);

            // Prepare contrast adjustment
            tileConfig = outputDirFFC_TP_Path+"/"+"TileConfiguration.registered.txt";
            std::string tileConfig2 = outputDirFFC_TP_Path + "/" + "TileConfiguration.registered2.txt";
            
            // Replace black tiles
            RemoveBlackTilesFromTileConfiguration(tileConfig, tileConfig2, blackTiles);
            //            std::string tileReplacement = "cp "+tileConfig2+" "+tileConfig;
            //            system(qPrintable(tileReplacement));
            
            
            std::string minInt = std::to_string(inputParamter.GetMinIntensity());
            std::string maxInt = std::to_string(inputParamter.GetMaxIntensity());
            
            std::string caProgram;
#ifdef __APPLE__
            caProgram = inputParamter.GetFiji() + " -batch " +
            inputParamter.GetScriptLocation() + "/RunContrastAdjustment_MacOsX.bsh" +
            " -in" + outputDirFFC_TP_Path + "/" +
            "=-out" + outputDirAC_TP_Path + "/" +
            "=-tile" + tileConfig2 +
            "=-min" + minInt +
            "=-max" + maxInt +
            "=-l1" + std::to_string(inputParamter.GetLambda1())+
            "=-l2" + std::to_string(inputParamter.GetLambda2());
            
#endif
            
#ifdef __linux__
            caProgram = "xvfb-run -a " + inputParamter.GetFiji() +
            " -Din=" + outputDirFFC_TP_Path + "/" +
            " -Dout=" + outputDirAC_TP_Path + "/" +
            " -Dtile=" + tileConfig2 +
            " -Dmin=" + minInt +
            " -Dmax=" + maxInt +
            " -Dl1=" + std::to_string(inputParamter.GetLambda1())+
            " -Dl2=" + std::to_string(inputParamter.GetLambda2()) +
            " -- --no-splash " + inputParamter.GetScriptLocation() + "/RunContrastAdjustment.bsh";
#endif
            system (&caProgram[0]);
            
            // Copy black tiles
            for (auto blackTile : blackTiles) {
              std::string copyCall = "cp " + outputDirFFC_TP_Path+ "/" + blackTile + " " + outputDirAC_TP_Path+ "/" + blackTile;
              system(&copyCall[0]);
            }
            
            // Actual Stitching
            
            std::string tileConfigCopy = outputDirAC_TP_Path+"/"+"TileConfiguration.txt";
            system(&("cp " + tileConfig2 + " " + tileConfigCopy)[0]);
            
#ifdef __APPLE__
            stitchProgram = inputParamter.GetFiji() + " -batch" +
            " " + inputParamter.GetScriptLocation() + "/gridStitching_UseFixedTileConfig_MacOsX.bsh"
            " -in" + outputDirAC_TP_Path +
            "=-out" + outputDirs.stitchedPath+ "/" +
            "=-tile" + "TileConfiguration.txt";
#endif
            
#ifdef __linux__
            stitchProgram = "xvfb-run -a " + inputParamter.GetFiji() +
            " -Din=" + outputDirAC_TP_Path +
            " -Dout=" + outputDirs.stitchedPath+ "/" +
            " -Dtile=" + "TileConfiguration.txt" +
            " -- --no-splash " + inputParamter.GetScriptLocation() + "/gridStitching_UseFixedTileConfig.bsh";
#endif
            system(&stitchProgram[0]);
            
            system(&("mv " + outputDirs.stitchedPath+"/"+"img_t1_z1_c1" + " " + outputDirs.stitchedPath+"/"+"t"+iPadded+"_preprocessed.tif")[0]);
            
            if (inputParamter.GetIntermediateFileDeletion()) {
              system(&("rm -rf "+ outputDirAC_TP_Path)[0]);
              system(&("rm -rf "+ outputDirFFC_TP_Path)[0]);
            }
            
          } else {
            std::cout << "No folder for time point " << i << " ----> Skip time point" << std::endl;
          }
        }
        
      }
      
      // Cleaning
      if (inputParamter.GetFlatField() != "") {
        Free_Array(flatField);
      }
      if (inputParamter.GetIntermediateFileDeletion()) {
        system(&("rm -rf "+ outputDirs.flatFieldPath)[0]);
        system(&("rm -rf "+ outputDirs.contrastAdjPath)[0]);
      }
      
    }
  }
  
}