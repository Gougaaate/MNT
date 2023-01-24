# MNT Project
This project generates a PPM image from a list of GPS coordinates and their corresponding depths. The image is created by first projecting the coordinates onto a Cartesian plane using the proj.h library. The coordinates are then used to create pixels, and hillshading is applied to give the image a 3D effect. The resulting image is then saved as a PPM file.

## Usage
To use this project, you will need to provide a list of GPS coordinates and their corresponding depths in a .txt format. The project can then be run from the command line, with the path to the input file as an argument.

## Copy code
Execute the build.sh file, or use this command :
    ```
    $ mkdir build; cd build; cmake ..; make```
    
The resulting PPM image will be saved in the same directory as the input file, with the 'MNT.ppm'.
To generate it : 

   
    $ ./create_raster *nameOfYourTxtFile* *pixel width*
    
Note : if the file name is "data.txt", type only :
```console
$ ./create_raster data 1000
```
  

## Input Format
The input file should be a plain text file, with one coordinate-depth pair per line. The coordinates should be in latitude and longitude, and the depth should be in meters. The values should be separated by a space.

## Source file required
A .txt file containing values separated by a blank :
latitude longitude depth   (GPS coordinates Lambert93)
latitude longitude depth
...
 ## Output Format
The output image will be a PPM file, which can be opened with most image viewers. The colors of the image will represent the depths of the coordinates, with darker colors representing deeper depths.

## Dependencies
This project requires the proj.h library for coordinate projection, the delaunator.hpp header to triangulate the image, and the C++ standard library for input/output and image generation.
You also need cmake to compile the project with the code line. Else, use g++.

