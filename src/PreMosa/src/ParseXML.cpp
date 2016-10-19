//
//  ParseXML.cpp
//  CMakeProject
//
//  Created by blasse on 7/25/13.
//
//

#include "ParseXML.h"

#include <iostream>
#include <QFile>

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
  XmlFile::XmlFile (QString xmlFileName) {
        
    QFile file(xmlFileName);
    if (!file.open(QIODevice::ReadOnly)) {
      std::cout << "Error opening XML: " << qPrintable(file.fileName()) << std::endl;

    }
    
    QDomDocument xmlDocument("parameterXML");
    if (!xmlDocument.setContent(&file)) {
      file.close();
      std::cout << "Error reading XML: " << qPrintable(file.fileName()) << std::endl;
    }
    file.close();
    
    QDomElement xmlElements = xmlDocument.documentElement();
    complete_ = true;
      
    ParseRequieredParameter(xmlElements);
    ParseNonRequieredParameter(xmlElements);
  }
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Destructor
   */
  XmlFile::~XmlFile () {
  }
  
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Parsing function of required parameters
   -> Individual elements are parsed in order to set the private parameters
   
   @param xmlElements:  all elements in the xml file
   */
  
  void XmlFile::ParseRequieredParameter(QDomElement xmlElements){
    
    // Read input folder
    QDomNodeList parameterNode = xmlElements.elementsByTagName("inputFolder");
    if (parameterNode.count() == 0 or parameterNode.at(0).toElement().text().length() == 0) {
      std::cout << "Error in XML: No input folder defined!" << std::endl;
      complete_ = false;
    } else {
      inputDir_ = parameterNode.at(0).toElement().text();
    }
    
    // Read output folder
    parameterNode = xmlElements.elementsByTagName("outputFolder");
    if (parameterNode.count() == 0 or parameterNode.at(0).toElement().text().length() == 0) {
      std::cout << "Error in XML: No output folder defined!" << std::endl;
      complete_ = false;
    } else {
      outputDir_ = parameterNode.at(0).toElement().text();
    }
    
    // Read folder name expression
    parameterNode = xmlElements.elementsByTagName("folderNameExpression");
    if (parameterNode.count() == 0 or parameterNode.at(0).toElement().text().length() == 0) {
      std::cout << "Error in XML: No folder name expression defined!" << std::endl;
      complete_ = false;
    } else {
      folderNameExpression_ = parameterNode.at(0).toElement().text();
    }
    
    // Read file name expression
    parameterNode = xmlElements.elementsByTagName("fileNameExpression");
    if (parameterNode.count() == 0 or parameterNode.at(0).toElement().text().length() == 0) {
      std::cout << "Error in XML: No file name expression defined!" << std::endl;
      complete_ = false;
    } else {
      fileNameExpression_ = parameterNode.at(0).toElement().text();
    }

    // Read path to the fiji scripts
    parameterNode = xmlElements.elementsByTagName("scriptLocation");
    if (parameterNode.count() == 0 or parameterNode.at(0).toElement().text().length() == 0) {
      std::cout << "Error in XML: No path to the fiji scripts defined!" << std::endl;
      complete_ = false;
    } else {
      scriptLocation_ = parameterNode.at(0).toElement().text();
    }
      
    // Read path to the fiji program
    parameterNode = xmlElements.elementsByTagName("fiji");
    if (parameterNode.count() == 0 or parameterNode.at(0).toElement().text().length() == 0) {
      std::cout << "Error in XML: No path to the fiji binary defined!" << std::endl;
      complete_ = false;
    } else {
      fiji_ = parameterNode.at(0).toElement().text();
    }
    
    // Parse time range
    parameterNode = xmlElements.elementsByTagName("tStart");
    if (parameterNode.count() == 0 or parameterNode.at(0).toElement().text().length() == 0) {
      std::cout << "Error in XML: No start time defined!" << std::endl;
      complete_ = false;
    }
    tStart_ = parameterNode.at(0).toElement().text().toInt();
    
    parameterNode = xmlElements.elementsByTagName("tEnd");
    if (parameterNode.count() == 0 or parameterNode.at(0).toElement().text().length() == 0) {
      std::cout << "Error in XML: No end time defined!" << std::endl;
      complete_ = false;
    }
    tEnd_ = parameterNode.at(0).toElement().text().toInt();
    
    
    // Parse grid size
    parameterNode = xmlElements.elementsByTagName("gridSizeX");
    if (parameterNode.count() == 0 or parameterNode.at(0).toElement().text().length() == 0) {
      std::cout << "Error in XML: No grid size in x direction defined!" << std::endl;
      complete_ = false;
    }
    gridSizeX_ = parameterNode.at(0).toElement().text().toInt();
    
    parameterNode = xmlElements.elementsByTagName("gridSizeY");
    if (parameterNode.count() == 0 or parameterNode.at(0).toElement().text().length() == 0) {
      std::cout << "Error in XML: No grid size in y direction defined!" << std::endl;
      complete_ = false;
    }
    gridSizeY_ = parameterNode.at(0).toElement().text().toInt();
    
  }
  
  
  
  ////////////////////////////////////////////////////////////////////////////////
  /*
   Parsing function of non-required parameters
   -> Individual elements are parsed in order to set the private parameters
      If an element is not specified in the xml file, a default value is assigned
   
   @param xmlElements:  all elements in the xml file
   */
  
  void XmlFile::ParseNonRequieredParameter(QDomElement xmlElements){
    
    // Parse verbose value
    QDomNodeList parameterNode = xmlElements.elementsByTagName("verbose");
    verbose_ = false;
    if (parameterNode.count() != 0) {
      if (parameterNode.at(0).toElement().text() == "true") {
        verbose_ = true;
      }
    }
    
    // Parse print height map value
    parameterNode = xmlElements.elementsByTagName("printHeightMap");
    printHeightMap_ = false;
    if (parameterNode.count() != 0) {
      if (parameterNode.at(0).toElement().text() == "true") {
      printHeightMap_ = true;
      }
    }
    
    // Parse intermediateFileDeletion value
    parameterNode = xmlElements.elementsByTagName("intermediateFileDeletion");
    intermediateFileDeletion_ = false;
    if (parameterNode.count() != 0) {
      if (parameterNode.at(0).toElement().text() == "true") {
      intermediateFileDeletion_ = true;
      }
    }
    
    // Parse radius value
    parameterNode = xmlElements.elementsByTagName("radius");
    if (parameterNode.count() != 0) {
      radius_ = parameterNode.at(0).toElement().text().toInt();
    } else {
      radius_ = 30;
    }
    
    // Parse distance value
    parameterNode = xmlElements.elementsByTagName("distance");
    if (parameterNode.count() != 0) {
      distance_ = parameterNode.at(0).toElement().text().toInt();
    } else {
      distance_ = 1;
    }
    
    // Parse threshold value
    parameterNode = xmlElements.elementsByTagName("threshold");
    if (parameterNode.count() != 0) {
      threshold_ = parameterNode.at(0).toElement().text().toInt();
    } else {
      threshold_ = 50;
    }
    
    // Parse threshold value
    parameterNode = xmlElements.elementsByTagName("overlap");
    if (parameterNode.count() != 0) {
      overlap_ = parameterNode.at(0).toElement().text().toInt();
    } else {
      overlap_ = 8;
    }
    
    // Read flat field file
    parameterNode = xmlElements.elementsByTagName("flatField");
    if (parameterNode.count() != 0) {
      flatField_ = parameterNode.at(0).toElement().text();
    }
    
    // Parse value for the first index of all tiles
    parameterNode = xmlElements.elementsByTagName("firstIndex");
    if (parameterNode.count() != 0) {
      firstIndex_ = parameterNode.at(0).toElement().text().toInt();
    } else {
        firstIndex_ = 0;
    }
    
    // Parse value for the minimal intensity being included into the contrast adjustment
    parameterNode = xmlElements.elementsByTagName("minIntensity");
    if (parameterNode.count() != 0) {
      minIntensity_ = parameterNode.at(0).toElement().text().toInt();
    } else {
      minIntensity_ = 0;
    }
    
    // Parse value for the maximal intensity being included into the contrast adjustment
    parameterNode = xmlElements.elementsByTagName("maxIntensity");
    if (parameterNode.count() != 0) {
      maxIntensity_ = parameterNode.at(0).toElement().text().toInt();
    } else {
      maxIntensity_ = 65535;
    }
    
    // Parse value for the number of pixel pairs sampled during the contrast adjustment
    parameterNode = xmlElements.elementsByTagName("numberOfPixelPairs");
    if (parameterNode.count() != 0) {
      numberOfSampledPairs_ = parameterNode.at(0).toElement().text().toInt();
    } else {
      numberOfSampledPairs_ = 100;
    }
    
    // Parse file path of the master tile configuration
    parameterNode = xmlElements.elementsByTagName("masterTileConfig");
    if (parameterNode.count() != 0) {
      masterTileConfig_ = parameterNode.at(0).toElement().text();
    }

    // Parse value for lambda1
    parameterNode = xmlElements.elementsByTagName("lambda1");
    if (parameterNode.count() != 0) {
      lamdba1_ = parameterNode.at(0).toElement().text().toDouble();
    } else {
      lamdba1_ = 0.1;
    }
    
    // Parse value for lambda2
    parameterNode = xmlElements.elementsByTagName("lambda2");
    if (parameterNode.count() != 0) {
      lamdba2_ = parameterNode.at(0).toElement().text().toDouble();
    } else {
      lamdba2_ = 0.1;
    }

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
  
  QString XmlFile::GetInputDir () const {
    return inputDir_;
  }
  
  QString XmlFile::GetOutputDir () const {
    return outputDir_;
  }
  
  QString XmlFile::GetFolderNameExpression () const {
    return folderNameExpression_;
  }
  
  QString XmlFile::GetFileNameExpression () const {
    return fileNameExpression_;
  }
    
  QString XmlFile::GetFiji () const {
    return fiji_;
  }

  QString XmlFile::GetScriptLocation () const {
    return scriptLocation_;
  }
  
  QString XmlFile::GetMasterTileConfig () const {
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
  
  QString XmlFile::GetFlatField () const {
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