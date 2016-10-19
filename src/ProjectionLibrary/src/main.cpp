//
//  main.cpp
//
//  Created by Corinna Blasse on 25.07.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include <string>
#include <time.h>


extern "C" {
  #include "utilities.h"
  #include "image.h"
  #include "mylib.h"
}

#include "ImageClass.h"
#include "Parameter.h"
#include "ProjectionAlgorithm.h"
#include "ShowArray.h"





#define DEVELOP


using namespace ImageClass;
using namespace Parameter;


/****************************************************************************************
 *                                                                                      *
 *  MAIN LEVEL                                                                          *
 *                                                                                      *
 ****************************************************************************************/

int main(int argc, char * argv[])
{
  time_t start,end;
  time (&start);

  static char *Spec[] = { "[!vg] [-r<int(20)>] [-t<int(50)>] [-d<int(0)>] [-l<int(1)>] [!h] [!s] [-acc] [-mi] [--version] [--range] <in:file> <out:file>", NULL };
  
  Process_Arguments(argc,argv,Spec,1);

  Parameter::ParaSet imgPara;
  imgPara.verbose    = Get_Int_Arg("-v");     //  Print out steps as they are completed
  imgPara.grid       = Get_Int_Arg("-g");        //  Paint grid on residue and marked images
  imgPara.radius     = Get_Int_Arg("-r");      //  The radius (-r) parameter value
  imgPara.accurate   = Get_Int_Arg("-acc");
  imgPara.maxInterpolation = Get_Int_Arg("-mi");
  imgPara.range = Get_Int_Arg("--range");
  imgPara.threshold  = Get_Int_Arg("-t");   //  The threshold (-s) parameter value
  imgPara.distance   = Get_Int_Arg("-d");    //  Distance parameter to define the reduced stack
  imgPara.printHeightMap   = Get_Int_Arg("-h");    //  Distance parameter to define the reduced stack
  imgPara.printRealHeightMap   = Get_Int_Arg("-s");    //  Distance parameter to define the reduced stack
  imgPara.version   = Get_Int_Arg("--version");
  imgPara.outputFolder = std::string(Get_String_Arg("out"));
  Image img = Image(Read_Image(Get_String_Arg("in"),0));
  imgPara.layer = Get_Int_Arg("-l");
  if (imgPara.layer < 1) {
    imgPara.layer = 1;
  }
  
  // Since the quadtree principle is used, the radius has to be even.
  if (imgPara.radius % 2 != 0) {
    imgPara.radius +=1;
  }
  img.SetGrid(imgPara.radius);
  
  if (imgPara.verbose) {
    std::cout << "Read in " << img.GetWidth() << "x" << img.GetHeight() << "x" << img.GetDepth() << " image" << std::endl;
    std::cout << "Image type: " << img.GetType() << " (0: uint8, 1:uint16)" << std::endl;
    std::cout << "Use an " << imgPara.radius << "x" << imgPara.radius << " window" << std::endl;
    std::cout << "Looking at the top " << imgPara.threshold << " pixels" << std::endl;
  }
  
  if (imgPara.version) {
    std::cout << "Projection algorithm 1.0" << std::endl << std::endl;
  }
  
  if (img.GetType() == UINT8_TYPE) {
    ProjectionAlgorithm::Project_8bit(img, imgPara);
  } else if (img.GetType() == UINT16_TYPE && imgPara.accurate) {
    ProjectionAlgorithm::Fast_Project_16bit(img, imgPara);
  } else if (img.GetType() == UINT16_TYPE) {
    //img.ScaleTo8bit();
    //projection = ProjectionAlgorithm::Project_8bit(img, imgPara);
    ProjectionAlgorithm::Fast_Project_16bit(img, imgPara);
  } else {
    std::cout << "Error: Image type is not supported!" << std::endl;
    img.~Image();
    return 0;
  }
  
    
  img.~Image();
  
  if (imgPara.verbose) {
    std::cout << "Array usage: " << Array_Usage() << std::endl;
  }
  
  time (&end);
  double dif = difftime (end,start);
  std::cout << "Time: " << dif << "s" << std::endl;
  
  return 0;
}
