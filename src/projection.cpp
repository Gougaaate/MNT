#include <proj.h>

void projection(const std::vector<double> &x_map, const std::vector<double> &y_map, std::vector<double> &new_xy)  {
    // Converts from GPS Point to cartesian Point
    PJ *P = proj_create_crs_to_crs(PJ_DEFAULT_CTX,"+proj=longlat +datum=WGS84",
                                   "+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",NULL);

    int n = x_map.size();
    for (int k = 0; k < n; k++) {
        PJ_COORD gps_coord, cartesian_coord;
        gps_coord.lpzt.lam = x_map[k];          // Longitude
        gps_coord.lpzt.phi = y_map[k];	        // Latitude

        cartesian_coord = proj_trans(P, PJ_FWD, gps_coord);
        new_xy.push_back(cartesian_coord.xy.x);                // Adding new values to the new (x,y) vectors
        new_xy.push_back(cartesian_coord.xy.y);
    }
}