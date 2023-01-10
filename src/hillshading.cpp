#ifndef __HILLSHADING_HPP__
#define __HILLSHADING_HPP__

#include <iostream>
#include "hillshading.hpp"
#include <cmath>
#include "display.hpp"


/**
 * @brief Calculates the shading value for each pixel in an elevation map
 * @param pixelVector 2D vector of pixel elevation values
 * @param height Number of rows in the pixel vector
 * @param width Number of columns in the pixel vector
 * @param pixelHeight Height of each pixel, in the same unit as the elevation values
 * @param pixelWidth Width of each pixel, in the same unit as the elevation values
 * @param azimuthDeg Azimuth angle of the light source, in degrees
 * @param heightDeg Height angle of the light source, in degrees
 * @return 2D vector of shading values, with the same dimensions as the input vector
 *
 * The function calculates a shading value for each pixel in the input vector,
 * based on the slope and aspect of the terrain at that pixel. The shading values
 * are calculated using the following formula:
 *
 *    shading = 255 * (cos(zenith) * cos(slope) + sin(zenith) * sin(slope) * cos(azimuth - aspect))
 *
 * where zenith, azimuth, slope, and aspect are the angles defined in the parameter description.
 * The function also call displayProgress function periodically to show the status of the computation.
 */

std::vector<std::vector<double>> hillshading(const std::vector<std::vector<double>> &pixelVector, int height, int width,  double pixelHeight, double pixelWidth, double azimuthDeg, double heightDeg) {
    const double zenithDeg = 90 - heightDeg;
    const double zenithRad = zenithDeg * M_PI / 180;                // Radian conversion
    double azimuthRad = 360 - azimuthDeg + 90;
    if (azimuthRad >= 360)
        azimuthRad = azimuthRad - 360;
    azimuthRad = azimuthRad * M_PI / 180;

    std::vector<std::vector<double>> pixelShading(height, std::vector<double>(width, 0));

    int count = 0;

    for (int i = 1; i < height -1; i++) {
        for (int j = 1; j < width -1; j++) {
            double a = pixelVector[i-1][j-1];            // e = center pixel, other letters = neighbours
            double b = pixelVector[i][j-1];
            double c = pixelVector[i+1][j-1];
            double d = pixelVector[i-1][j];
            double e = pixelVector[i][j];
            double f = pixelVector[i+1][j];
            double g = pixelVector[i-1][j+1];
            double h = pixelVector[i][j+1];
            double k = pixelVector[i+1][j+1];        // Not i because of the for (int i = ...)

            double dz_dx = ((c+k+(2*f))-(a+g+(2*d)))/(8 * pixelWidth);         // Derivatives
            double dz_dy = ((g+k+(2*h))-(a+c+(2*b)))/(8 * pixelHeight);

            double slopeRad = atan(e*sqrt((dz_dx*dz_dx)+(dz_dy*dz_dy)));
            double aspect_rad;

            if (dz_dx!=0)
            {
                aspect_rad = atan2(dz_dy,-1*dz_dx);
                if(aspect_rad<0)                                        // Modulo 2pi to avoid negative values
                    aspect_rad = 2*M_PI + aspect_rad;
            }
            else                // Indexing conditions (from the "how does Hillshading work" website)
            {
                if (dz_dy>0)
                    aspect_rad = M_PI/2;
                else
                {
                    if (dz_dy<0)
                        aspect_rad = 2*M_PI - (M_PI/2);
                    else
                        aspect_rad = 0;
                }
            }
            pixelShading[i][j] = 255.0 * ((cos(zenithRad) * cos(slopeRad)) + (sin(zenithRad) * sin(slopeRad) * cos(azimuthRad - aspect_rad)));      // Shading factor

            if (pixelShading[i][j]<0)           // Indexing
                pixelShading[i][j] = 0;
            if (pixelShading[i][j]>255)
                pixelShading[i][j] = 255;
            count ++;
            if (count % 500 == 0)                                       // Displaying the status
                displayProgress(100*count/((height-1)*(width-1)));
        }
    }
    return pixelShading;
}

#endif