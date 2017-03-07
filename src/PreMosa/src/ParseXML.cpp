//
//  ParseXML.cpp
//  CMakeProject
//
//  Created by blasse on 7/25/13.
//
//

#include "ParseXML.h"

#include <fstream>
#include <iostream>
#include <vector>
using namespace std;

namespace PreprocessingPipeline {
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Constructor (default)
   */
  
  XmlFile::XmlFile () {
  }
  
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Constructor with xml file path as input
   -> Parameter are parsed automatically
   
   @param xmlFileName:  file path to the xml file
   */
  XmlFile::XmlFile (std::string xmlFileName) {
    
    TiXmlDocument doc( const_cast<char *>(&xmlFileName[0]) );
    bool loadOkay = doc.LoadFile();
    if (!loadOkay) {
      std::cout << "Error reading XML: " << xmlFileName << std::endl;
    }
    TiXmlHandle xmlElements(&doc);
    complete_ = true;
    SetDefaultParameters();
    ParseRequieredParameter(doc);
  }
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Destructor
   */
  XmlFile::~XmlFile () {
  }
  
  
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Method to verify if all required parameters are set
   */
  void XmlFile::SetDefaultParameters () {
    gridSizeX_ = -1;
    gridSizeY_ = -1;
    tStart_ = -1;
    tEnd_ = -1;
    verbose_ = false;
    radius_ = 30;
    distance_ = 1;
    lamdba1_ = 0.1;
    lamdba2_ = 0.1;
    printHeightMap_ = false;
    intermediateFileDeletion_ = false;
    threshold_ = 50;
    overlap_ = 8;
    firstIndex_ = 0;
    minIntensity_ = 0;
    maxIntensity_ = 65535;
    numberOfSampledPairs_ = 100;
  }
  
  
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
   Method to verify if all required parameters are set
   */
  void XmlFile::VerifyRequiredParameters (bool hasStitchedTag) {
    
    if (inputDir_ == "") {
      std::cout << "Error in XML: No input folder defined!" << std::endl;
      complete_ = false;
    }
    if (outputDir_ == "") {
      std::cout << "Error in XML: No output folder defined!" << std::endl;
      complete_ = false;
    }
    if (folderNameExpression_ == "") {
      std::cout << "Error in XML: No folder name expression defined!" << std::endl;
      complete_ = false;
    }
    if (fileNameExpression_ == "") {
      std::cout << "Error in XML: No file name expression defined!" << std::endl;
      complete_ = false;
    }
    if (tStart_ == -1) {
      std::cout << "Error in XML: No start time defined!" << std::endl;
      complete_ = false;
    }
    if (tEnd_ == -1) {
      std::cout << "Error in XML: No end time defined!" << std::endl;
      complete_ = false;
    }
    if (fiji_ == "") {
      std::cout << "Error in XML: No path to the fiji binary defined!" << std::endl;
      complete_ = false;
    }
    if (scriptLocation_ == "") {
      std::cout << "Error in XML: No path to the fiji scripts defined!" << std::endl;
      complete_ = false;
    }
    if (gridSizeX_ == -1) {
      std::cout << "Error in XML: No grid size in x direction defined!" << std::endl;
      complete_ = false;
    }
    if (gridSizeY_ == -1) {
      std::cout << "Error in XML: No grid size in y direction defined!" << std::endl;
      complete_ = false;
    }
    if (hasStitchedTag == false) {
      std::cout << "Error in XML: No stitching element in the required element defined!" << std::endl;
      complete_ = false;
    }
    
    if (!FileExists(inputDir_)) {
      std::cout << "Error : Input folder does not exist!" << std::endl;
      complete_ = false;
    }
    if (!FileExists(scriptLocation_)) {
      std::cout << "Error : Fiji script location does not exist!" << std::endl;
      complete_ = false;
    }
    if (flatField_ != "") {
      if (!FileExists(flatField_)) {
        std::cout << "Error : Flat field image does not exist!" << std::endl;
        complete_ = false;
      }
    }
  }
  
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Parsing function of required parameters
   -> Individual elements are parsed in order to set the private parameters
   
   @param xmlElements:  all elements in the xml file
   */
  void XmlFile::ParseRequieredParameter(TiXmlDocument &xmlDoc){
    
    TiXmlElement * premosaTag = xmlDoc.FirstChildElement("premosa");
    if (!premosaTag) {
      std::cout << "Error in XML: The strutcure is not right (<premosa> element is missing). " << std::endl;
      complete_ = false;
      return;
    }
    
    bool hasStitchingTag = false;
    
    for ( TiXmlElement* parameter = premosaTag->FirstChildElement( );
         parameter;
         parameter = parameter->NextSiblingElement( ) ) {
      
      // Parsing of the parameters
      if (strcmp(parameter->Value( ), "required") == 0) {
        
        for ( TiXmlElement* requiredParameter = parameter->FirstChildElement( );
             requiredParameter;
             requiredParameter = requiredParameter->NextSiblingElement( ) ) {
          
          if (strcmp(requiredParameter->Value( ), "inputFolder") == 0) {
            if (requiredParameter->GetText() != NULL) {
              inputDir_ = requiredParameter->GetText();
            }
          } else if (strcmp(requiredParameter->Value( ), "outputFolder") == 0) {
            if (requiredParameter->GetText() != NULL) {
              outputDir_ = requiredParameter->GetText();
            }
          } else if (strcmp(requiredParameter->Value( ), "folderNameExpression") == 0) {
            if (requiredParameter->GetText() != NULL) {
              folderNameExpression_ = requiredParameter->GetText();
            }
          } else if (strcmp(requiredParameter->Value( ), "fileNameExpression") == 0) {
            if (requiredParameter->GetText() != NULL) {
              fileNameExpression_ = requiredParameter->GetText();
            }
          } else if (strcmp(requiredParameter->Value( ), "tStart") == 0) {
            if (requiredParameter->GetText() != NULL) {
              tStart_ = atoi(requiredParameter->GetText());
            }
          } else if (strcmp(requiredParameter->Value( ), "tEnd") == 0) {
            if (requiredParameter->GetText() != NULL) {
              tEnd_ = atoi(requiredParameter->GetText());
            }
          } else if (strcmp(requiredParameter->Value( ), "fiji") == 0) {
            if (requiredParameter->GetText() != NULL) {
              fiji_ = requiredParameter->GetText();
            }
          } else if (strcmp(requiredParameter->Value( ), "scriptLocation") == 0) {
            if (requiredParameter->GetText() != NULL) {
              scriptLocation_ = requiredParameter->GetText();
            }
          } else if (strcmp(requiredParameter->Value( ), "stitching") == 0) {
            hasStitchingTag = true;
            
            for ( TiXmlElement* stitchingParameter = requiredParameter->FirstChildElement( );
                 stitchingParameter;
                 stitchingParameter = stitchingParameter->NextSiblingElement( ) ) {
              
              if (strcmp(stitchingParameter->Value( ), "gridSizeX") == 0) {
                if (stitchingParameter->GetText() != NULL) {
                  gridSizeX_ = atoi(stitchingParameter->GetText());
                }
              } else if (strcmp(stitchingParameter->Value( ), "gridSizeY") == 0) {
                if (stitchingParameter->GetText() != NULL) {
                  gridSizeY_ = atoi(stitchingParameter->GetText());
                }
              } else if (strcmp(stitchingParameter->Value( ), "firstIndex") == 0) {
                if (stitchingParameter->GetText() != NULL) {
                  firstIndex_ = atoi(stitchingParameter->GetText());
                }
              }
            }
          }
        }
        
      } else if (strcmp(parameter->Value( ), "optional")== 0) {
        
        for ( TiXmlElement* optionalParameter = parameter->FirstChildElement( );
             optionalParameter;
             optionalParameter = optionalParameter->NextSiblingElement( ) ) {
          
          if (strcmp(optionalParameter->Value( ), "verbose") == 0) {
            if (optionalParameter->GetText() != NULL) {
              if (strcmp(optionalParameter->GetText(), "true") == 0) {
                verbose_ = true;
              }
            }
          } else if (strcmp(optionalParameter->Value( ), "projection") == 0) {
            
            for ( TiXmlElement* projectionParameter = optionalParameter->FirstChildElement( );
                 projectionParameter;
                 projectionParameter = projectionParameter->NextSiblingElement( ) ) {
              
              if (strcmp(projectionParameter->Value( ), "radius") == 0) {
                if (projectionParameter->GetText() != NULL) {
                  radius_ = atoi(projectionParameter->GetText());
                }
              } else if (strcmp(projectionParameter->Value( ), "distance") == 0) {
                if (projectionParameter->GetText() != NULL) {
                  distance_ = atoi(projectionParameter->GetText());
                }
              } else if (strcmp(projectionParameter->Value( ), "threshold") == 0) {
                if (projectionParameter->GetText() != NULL) {
                  threshold_ = atoi(projectionParameter->GetText());
                }
              } else if (strcmp(projectionParameter->Value( ), "printHeightMap") == 0) {
                if (projectionParameter->GetText() != NULL) {
                  if (strcmp(projectionParameter->GetText(), "true") == 0) {
                    printHeightMap_ = true;
                  }
                }
              }
            }
          } else if (strcmp(optionalParameter->Value( ), "stitching") == 0) {
            
            for ( TiXmlElement* stitchingParameter = optionalParameter->FirstChildElement( );
                 stitchingParameter;
                 stitchingParameter = stitchingParameter->NextSiblingElement( ) ) {
              
              if (strcmp(stitchingParameter->Value( ), "overlap") == 0) {
                if (stitchingParameter->GetText() != NULL) {
                  overlap_ = atoi(stitchingParameter->GetText());
                }
              } else if (strcmp(stitchingParameter->Value( ), "masterTileConfig") == 0) {
                if (stitchingParameter->GetText() != NULL) {
                  masterTileConfig_ = stitchingParameter->GetText();
                }
              }
            }
          } else if (strcmp(optionalParameter->Value( ), "contrastOptimization") == 0) {
            
            for ( TiXmlElement* contrastParameter = optionalParameter->FirstChildElement( );
                 contrastParameter;
                 contrastParameter = contrastParameter->NextSiblingElement( ) ) {
              
              if (strcmp(contrastParameter->Value( ), "flatField") == 0) {
                if (contrastParameter->GetText() != NULL) {
                  flatField_ = contrastParameter->GetText();
                }
              } else if (strcmp(contrastParameter->Value( ), "minIntensity") == 0) {
                if (contrastParameter->GetText() != NULL) {
                  minIntensity_ = atoi(contrastParameter->GetText());
                }
              } else if (strcmp(contrastParameter->Value( ), "maxIntensity") == 0) {
                if (contrastParameter->GetText() != NULL) {
                  maxIntensity_ = atoi(contrastParameter->GetText());
                }
              } else if (strcmp(contrastParameter->Value( ), "numberOfPixelPairs") == 0) {
                if (contrastParameter->GetText() != NULL) {
                  numberOfSampledPairs_ = atoi(contrastParameter->GetText());
                }
              } else if (strcmp(contrastParameter->Value( ), "lambda1") == 0) {
                if (contrastParameter->GetText() != NULL) {
                  lamdba1_ = atof(contrastParameter->GetText());
                }
              } else if (strcmp(contrastParameter->Value( ), "lambda2") == 0) {
                if (contrastParameter->GetText() != NULL) {
                  lamdba2_ = atof(contrastParameter->GetText());
                }
              }
            }
          } else if (strcmp(optionalParameter->Value( ), "intermediateFileDeletion") == 0) {
            if (optionalParameter->GetText() != NULL) {
              if (strcmp(optionalParameter->GetText(), "true") == 0) {
                intermediateFileDeletion_ = true;
              }
            }
          }
        }
        
      } else {
        std::cout << "Error in XML: The strutcure is not right (<required> and/ or <optional> element is missing). " << std::endl;
        complete_ = false;
      }
    }
    
    // Verify if all required parameter are added
    VerifyRequiredParameters(hasStitchingTag);
  }
  
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Return functions
   */
  
  bool XmlFile::GetComplete () const {
    return complete_;
  }
  
  bool XmlFile::GetPrintHeightMap () const {
    return printHeightMap_;
  }
  
  bool XmlFile::GetIntermediateFileDeletion () const {
    return intermediateFileDeletion_;
  }
  
  std::string XmlFile::GetInputDir () const {
    return inputDir_;
  }
  
  std::string XmlFile::GetOutputDir () const {
    return outputDir_;
  }
  
  std::string XmlFile::GetFolderNameExpression () const {
    return folderNameExpression_;
  }
  
  std::string XmlFile::GetFileNameExpression () const {
    return fileNameExpression_;
  }
  
  std::string XmlFile::GetFiji () const {
    return fiji_;
  }
  
  std::string XmlFile::GetScriptLocation () const {
    return scriptLocation_;
  }
  
  std::string XmlFile::GetMasterTileConfig () const {
    return masterTileConfig_;
  }
  
  int XmlFile::GetTStart () const {
    return tStart_;
  }
  
  int XmlFile::GetTEnd () const {
    return tEnd_;
  }
  
  int XmlFile::GetGridSizeX () const {
    return gridSizeX_;
  }
  
  int XmlFile::GetGridSizeY () const {
    return gridSizeY_;
  }
  
  int XmlFile::GetFirstIndex () const {
    return firstIndex_;
  }
  
  int XmlFile::GetOverlap() const {
    return overlap_;
  }
  
  int XmlFile::GetVerbose () const {
    return verbose_;
  }
  
  std::string XmlFile::GetFlatField () const {
    return flatField_;
  }
  
  int XmlFile::GetRadius () const {
    return radius_;
  }
  
  int XmlFile::GetDistance () const {
    return distance_;
  }
  
  int XmlFile::GetThreshold () const {
    return threshold_;
  }
  
  int XmlFile::GetMinIntensity () const {
    return minIntensity_;
  }
  
  int XmlFile::GetMaxIntensity () const {
    return maxIntensity_;
  }
  
  int XmlFile::GetNumberOfSampledPairs () const {
    return numberOfSampledPairs_;
  }
  
  double XmlFile::GetLambda1 () const {
    return lamdba1_;
  }
  
  double XmlFile::GetLambda2 () const {
    return lamdba2_;
  }
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Functions to afterwards change the minIntensity and maxIntensity
   */
  
  void XmlFile::SetMinIntensity(int minIntensity) {
    minIntensity_ = minIntensity;
  }
  
  void XmlFile::SetMaxIntensity(int maxIntensity) {
    maxIntensity_ = maxIntensity;
  }
  
}