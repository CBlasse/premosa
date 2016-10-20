PreMosa
------------------------
The imaging of large tissues or specimen with fluorescent microscopy typically yields large 3D image mosaics. In many experimental setups, however, the used markers localize only along a 2D manifold (surface) within the tissue. To simplify downstream analyses, one can therefore reduce the image data to this 2D manifold. PreMosa is a fully automatic pipeline to perform this task and three additional tasks that are also essential in the preprocessing of 3D image mosaics:

1) Extraction and projection of the 2D surface
2) Correction of uneven illumination artifacts
3) Stitching of the 2D mosaic planes
4) Adjusting contrast and brightness of each 2D plane

For each image mosaic, the preprocessing results a single, large 2D image that presents the stained surface of interest without illumination artifacts.

Authors and Contributors
-------------------------

The PREMOSA project was initiated by Corinna Blasse and Gene Myers. Corinna Blasse wrote the algorithms for extracting 2D surfaces in 3D volumes, for correcting the illumination artifacts and for combining the different tasks into one pipeline. Stephan Saalfeld wrote the Fiji plugin to adjust the contrast and brightness of the 2D mosaic planes.

Installation
-------------------------

The pipeline has been only developed for Linux/ Mac OSX.

Prerequisites:
  - Linux / Mac OSX system
  - Cmake (get it here)
  - Qt (version 5.6, get it here)
  - Fiji (get it here)

Installation:

1) Download the source code by using either the download page or the git repository: git clone https://github.com/CBlasse/premosa.git

2) Go to the project folder and run: ./Build.sh -p Qt_prefix_path_to_cmake_configurations
  This step will generate executable programs in the ./bin/ folder. It is thereby necessary to add the Qt prefix path for finding the cmake configuration files (e.g. /usr/local/Qt5.6/5.6/clang_64/lib/cmake).

3) Install the contrast adjustment plugin in Fiji. To do so, copy the contrastAdjustment.jar into the plugins folder of your Fiji distribution.



Usage
-------------------------
Please have a look at https://cblasse.github.io/premosa/usage.html