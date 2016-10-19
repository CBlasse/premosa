//
//  Pipeline.cpp
//  CMakeProject
//
//  Created by blasse on 7/26/13.
//
//

#include "Pipeline.h"

#include <math.h>
#include <vector>

#include <QFile>
#include <QString>
#include <QSysInfo>
#include <QTextStream>

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
     Creation of a tile configuration using a master tile configuration
     
     @param inputParamter:    struct with parameters
     @param tileConfigPath:   file name of the new tile configuration file
     @param fileName:         name template of the input images
     */
    
    
    void CreateTileConfiguration (XmlFile inputParamter, QString tileConfigPath, QString fileName) {
      
      QFile masterFile(inputParamter.GetMasterTileConfig());
      
      if (!masterFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
        std::cout << "Master tile configuration file cannot be opened! A new configuration will be computed!" << std::endl;
        return;
      }
      
      int paddingNumber = inputParamter.GetFileNameExpression().indexOf("}") - inputParamter.GetFileNameExpression().indexOf("{") - 1;  // determination of the number of digits of the tile indication
      
      QTextStream in(&masterFile);
      QString lines = in.readAll();
      
      for (int i=inputParamter.GetFirstIndex(); i<=inputParamter.GetGridSizeX()*inputParamter.GetGridSizeY(); i++) {  // iteration through the number of tile
        
        QString expression = '{' + QString('i').rightJustified(paddingNumber, 'i')+'}'; // replace the current tile number
        QString tile = QString::number(i).rightJustified(paddingNumber, '0');
        QString tileReplacement = "%TILE" + QString::number(i) +"%";
        QString newFileName = fileName;
        newFileName.replace(expression, tile);
        
        lines = lines.replace(tileReplacement, newFileName);
      }
      
      QFile tileConfig(tileConfigPath); // writing of the new tile configuration file
      tileConfig.open(QIODevice::WriteOnly | QIODevice::Text);
      QTextStream dataTileConfig(&tileConfig);
      
      dataTileConfig << lines;
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
    
    void RemoveBlackTilesFromTileConfiguration (QString tileConfigPath,QString tileConfigPath2, std::vector<std::string> &blackTiles) {
      
      QFile masterFile(tileConfigPath);
      QFile masterFile2(tileConfigPath2);
      
      masterFile.open(QIODevice::ReadOnly | QIODevice::Text) ;
      masterFile2.open(QIODevice::WriteOnly | QIODevice::Text);
      QTextStream in(&masterFile);
      QTextStream out(&masterFile2);
      QString line;
      
      while (!in.atEnd())
        {
        line = in.readLine();
        std::string fileName = line.split(";")[0].toStdString();
        if (std::find(blackTiles.begin(), blackTiles.end(), fileName) == blackTiles.end()) {    // if the file name is not included in the black tile vector, then the line is added to the new tile configuration
          out << line << "\n";
        }
        }
      
      masterFile.close();
      masterFile2.close();
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
      
      QDir inputDir(inputParamter.GetInputDir());
      OutputDirs outputDirs = MakeOutputDirs (inputParamter.GetOutputDir());  // Creation of the required output folders
      
      int paddingNumber = inputParamter.GetFolderNameExpression().indexOf("}") - inputParamter.GetFolderNameExpression().indexOf("{") - 1;    // Parsing the number of digits used to encode the time point
      QString expression = '{' + QString('i').rightJustified(paddingNumber, 'i')+'}';
      
      
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
        flatField = Read_Image(const_cast<char * > (qPrintable(inputParamter.GetFlatField())),0);
        Range_Bundle flatFieldRange;
        Array_Range(&flatFieldRange, flatField);
        flatFieldMaximum = ceil(flatFieldRange.maxval.fval);
      }
      
      
      if (inputDir.exists()) {
        
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
          QString iPadded = QString::number(inputParamter.GetTStart()).rightJustified(paddingNumber, '0');
          QString folderName (inputParamter.GetFolderNameExpression());
          folderName.replace(expression, iPadded);
          QDir exampleDir (inputDir.absolutePath()+QDir::separator()+folderName);
          
          Array * exampleImage = Read_Image(const_cast<char * > (qPrintable(exampleDir.entryInfoList(QDir::Files)[0].filePath())),0);
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
          
          QString iPadded = QString::number(i).rightJustified(paddingNumber, '0');
          QString folderName (inputParamter.GetFolderNameExpression());
          folderName.replace(expression, iPadded);
          
          std::cout << "Processed time point: " << qPrintable(inputDir.absolutePath()+QDir::separator()+folderName) << std::endl;
          QDir timePointDir(inputDir.absolutePath()+QDir::separator()+folderName);
          
          QStringList filters;
          filters << "*.tif";
          timePointDir.setNameFilters(filters);
          std::vector<std::string> blackTiles;
          
          if (timePointDir.exists()) {
            QFileInfoList entries = timePointDir.entryInfoList(QDir::Files);
            
            // Create subdirectory
            QDir outputDirFFC_TP(outputDirs.flatFieldDir.absolutePath()+QDir::separator()+folderName);
            if (not outputDirFFC_TP.exists()) {
              outputDirFFC_TP.mkdir(outputDirs.flatFieldDir.absolutePath()+QDir::separator()+folderName);
            }
            
            foreach ( QFileInfo entryInfo, entries){
              
              // Read the image
              QString path = entryInfo.absoluteFilePath();

              Array * image = Read_Image(const_cast<char * > (qPrintable(entryInfo.filePath())),0);
              
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
                  projectionAlgorithm.DrawInterpolatedHeightMap(const_cast<char * > (qPrintable(outputDirFFC_TP.path()+QDir::separator()+entryInfo.fileName())));
                }
              }

              
              // Review if the image features significant signals and is not a 'black tile'
              CheckForBlackTile(projection, entryInfo.fileName().toStdString(), blackTiles);
              
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
                Write_Image(const_cast<char * > (qPrintable(outputDirFFC_TP.path()+QDir::separator()+entryInfo.fileName())), correctedImg, DONT_PRESS);
                Free_Array(correctedImg);
              } else {
                Write_Image(const_cast<char * > (qPrintable(outputDirFFC_TP.path()+QDir::separator()+entryInfo.fileName())), projection, DONT_PRESS);
              }
              
              // Cleaning (projection is freed by the destructor of the ProjectionAlgorithm)
              // Free_Array(projection);
              Free_Array(image);
            }
            
            // Folder for the contrast adjustment
            QDir outputDirAC_TP(outputDirs.constrastAdjDir.absolutePath()+QDir::separator()+folderName);
            if (not outputDirAC_TP.exists()) {
              outputDirAC_TP.mkdir(outputDirs.constrastAdjDir.absolutePath()+QDir::separator()+folderName);
            }
            
            QString fileName (inputParamter.GetFileNameExpression());
            fileName.replace("TIME", iPadded);
            
            
            // Mosaic stitching to determine the tile configuration
            QString tileConfig = outputDirFFC_TP.absolutePath() + QDir::separator() + "TileConfiguration.txt";
            if (inputParamter.GetMasterTileConfig() != "") {
              CreateTileConfiguration(inputParamter, tileConfig, fileName);
            }
            QFile tileFile (tileConfig);
            QString stitchProgram;
            
            // Fiji call (depends on the presence of a master tile configuration and the system(Linux/OsX) )
            if (!tileFile.exists()) {
              
#ifdef Q_WS_MAC
              stitchProgram = inputParamter.GetFiji() + " -batch" +
              " " + inputParamter.GetScriptLocation() + "/gridStitching_GetTileConfig_MacOsX.bsh" +
              " -i" + outputDirFFC_TP.absolutePath() +
              "=-x" + QString::number(inputParamter.GetGridSizeX()) +
              "=-y" + QString::number(inputParamter.GetGridSizeY()) +
              "=-l" +QString::number(inputParamter.GetOverlap())+
              "=-n" + fileName +
              "=-f"+ QString::number(inputParamter.GetFirstIndex());
#endif
              
#ifdef __linux__
              stitchProgram = "xvfb-run -a " + inputParamter.GetFiji() +
              " -Di=" + outputDirFFC_TP.absolutePath() +
              " -Dx=" + QString::number(inputParamter.GetGridSizeX()) +
              " -Dy=" + QString::number(inputParamter.GetGridSizeY()) +
              " -Dl=" +QString::number(inputParamter.GetOverlap())+
              " -Dn=" + fileName +
              " -Df="+ QString::number(inputParamter.GetFirstIndex()) +
              " -- --no-splash " + inputParamter.GetScriptLocation() + "/gridStitching_GetTileConfig.bsh";
#endif
              
            } else {
#ifdef Q_WS_MAC
              stitchProgram = inputParamter.GetFiji() + " -batch" +
              " " + inputParamter.GetScriptLocation() + "/gridStitching_UseTileConfig_MacOsX.bsh"+
              " -in" + outputDirFFC_TP.absolutePath() +
              "=-out" + outputDirFFC_TP.absolutePath()+QDir::separator();
              
#endif
              
#ifdef __linux__
              stitchProgram = "xvfb-run -a " + inputParamter.GetFiji() +
              " -Din=" + outputDirFFC_TP.absolutePath() +
              " -Dout=" + outputDirFFC_TP.absolutePath()+QDir::separator() +
              " -- --no-splash " + inputParamter.GetScriptLocation() + "/gridStitching_UseTileConfig.bsh";
#endif
            }
            system(qPrintable(stitchProgram));
            
            
            // Prepare contrast adjustment
            tileConfig = outputDirFFC_TP.absolutePath()+QDir::separator()+"TileConfiguration.registered.txt";
            QString tileConfig2 = outputDirFFC_TP.absolutePath() + QDir::separator() + "TileConfiguration.registered2.txt";
            
            // Replace black tiles
            RemoveBlackTilesFromTileConfiguration(tileConfig, tileConfig2, blackTiles);
            //            QString tileReplacement = "cp "+tileConfig2+" "+tileConfig;
            //            system(qPrintable(tileReplacement));
            
            
            QString minInt = QString::number(inputParamter.GetMinIntensity());
            QString maxInt = QString::number(inputParamter.GetMaxIntensity());
            
            QString caProgram;
#ifdef Q_WS_MAC
            //            caProgram = inputParamter.GetFiji() +
            //            " --run \"Contrast Adjustment\" \"folder=" + outputDirFFC_TP.absolutePath() + QDir::separator() +
            //            " output_folder=" + outputDirAC_TP.absolutePath() + QDir::separator() +
            //            " tile_configuration=" + tileConfig +
            //            " minimum_intensity=" + minInt +
            //            " maximum_intensity=" + maxInt +
            //            " number_of_samples=100 lambda_1=0.10 lambda_2=0.10\"";
            
            caProgram = inputParamter.GetFiji() + " -batch " +
            inputParamter.GetScriptLocation() + "/RunContrastAdjustment_MacOsX.bsh" +
            " -in" + outputDirFFC_TP.absolutePath() + QDir::separator() +
            "=-out" + outputDirAC_TP.absolutePath() + QDir::separator() +
            "=-tile" + tileConfig2 +
            "=-min" + minInt +
            "=-max" + maxInt;
            
#endif
            
#ifdef __linux__
            caProgram = "xvfb-run -a " + inputParamter.GetFiji() +
            " -Din=" + outputDirFFC_TP.absolutePath() + QDir::separator() +
            " -Dout=" + outputDirAC_TP.absolutePath() + QDir::separator() +
            " -Dtile=" + tileConfig2 +
            " -Dmin=" + minInt +
            " -Dmax=" + maxInt +
            " -- --no-splash " + inputParamter.GetScriptLocation() + "/RunContrastAdjustment.bsh";
#endif
            system(qPrintable(caProgram));
            
            // Copy black tiles
            for (auto blackTile : blackTiles) {
              QString copyCall = "cp " + outputDirFFC_TP.absolutePath()+ QDir::separator() + QString::fromStdString(blackTile) + " " + outputDirAC_TP.absolutePath()+ QDir::separator() + QString::fromStdString(blackTile);
              system(qPrintable(copyCall));
            }
            
            // Actual Stitching
            
            QString tileConfigCopy = outputDirAC_TP.absolutePath()+QDir::separator()+"TileConfiguration.txt";
            system(qPrintable("cp " + tileConfig2 + " " + tileConfigCopy));
            
#ifdef Q_WS_MAC
            stitchProgram = inputParamter.GetFiji() + " -batch" +
            " " + inputParamter.GetScriptLocation() + "/gridStitching_UseFixedTileConfig_MacOsX.bsh"
            " -in" + outputDirAC_TP.absolutePath() +
            "=-out" + outputDirs.stitchedDir.absolutePath()+QDir::separator() +
            "=-tile" + "TileConfiguration.txt";
#endif
            
#ifdef __linux__
            stitchProgram = "xvfb-run -a " + inputParamter.GetFiji() +
            " -Din=" + outputDirAC_TP.absolutePath() +
            " -Dout=" + outputDirs.stitchedDir.absolutePath()+QDir::separator() +
            " -Dtile=" + "TileConfiguration.txt" +
            " -- --no-splash " + inputParamter.GetScriptLocation() + "/gridStitching_UseFixedTileConfig.bsh";
#endif
            system(qPrintable(stitchProgram));
            
            system(qPrintable("mv " + outputDirs.stitchedDir.absolutePath()+QDir::separator()+"img_t1_z1_c1" + " " + outputDirs.stitchedDir.absolutePath()+QDir::separator()+"t"+iPadded+"_preprocessed.tif"));
            
            if (inputParamter.GetIntermediateFileDeletion()) {
              system(qPrintable("rm -rf "+ outputDirAC_TP.absolutePath()));
              system(qPrintable("rm -rf "+ outputDirFFC_TP.absolutePath()));
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
        system(qPrintable("rm -rf "+ outputDirs.flatFieldDir.absolutePath()));
        system(qPrintable("rm -rf "+ outputDirs.constrastAdjDir.absolutePath()));
      }
      
    }
  }
  
}