#ifndef __PROJECTION_H__
#define __PROJECTION_H__

#include <proj.h>
#include <vector>
#include "display.hpp"

/**
 * @brief Transforms GPS coordinates to cartesian coordinates
 * @param x_map  Vector of x (longitude) GPS coordinates
 * @param y_map  Vector of y (latitude) GPS coordinates
 * @param new_xy Vector of transformed cartesian coordinates (x,y)
 *
 * The function uses the PROJ library to convert GPS coordinates in WGS84 (longitude and latitude) to cartesian coordinates (x,y) in a Lambert Conformal Conic projection.
 * It takes two vectors of GPS coordinates, and return a vector of cartesian coordinates, Each GPS coordinates is transformed and added to the new_xy vector.
 */

void GPS_to_cartesian(const std::vector<double> &x_map, const std::vector<double> &y_map, std::vector<double> &new_xy)  {
    PJ *P = proj_create_crs_to_crs(PJ_DEFAULT_CTX,"+proj=longlat +datum=WGS84",
                                   "+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
                                   nullptr);

    int n = (int) x_map.size();
    for (int k = 0; k < n; k++) {
        PJ_COORD gps_coord, cartesian_coord;
        gps_coord.lpzt.lam = x_map[k];          // Longitude
        gps_coord.lpzt.phi = y_map[k];	        // Latitude

        cartesian_coord = proj_trans(P, PJ_FWD, gps_coord);
        new_xy.push_back(cartesian_coord.xy.x);                // Adding new values to the new (x,y) vectors
        new_xy.push_back(cartesian_coord.xy.y);
    }
}

#endif