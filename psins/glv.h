#pragma once
#include <Eigen/Dense>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace Eigen;

class GLV {
public:
    double deg, sec, min, hour, ppm;
    double dph, dpsh, arcmin, arcsec;
    double Re, wie, f, e2, ge, g0;
    double m, beta, beta1;
    double ug, Hz, ugpsHz, ugpg2, nm;
    double ts;
    GLV() {
        ts = 1.0 / 400.0;
        deg = M_PI / 180.0;
        sec = 1.0;
        min = 60.0;
        hour = 3600.0;
        ppm = 1e-6;
        
        dph = deg / hour;       
        dpsh = deg / sqrt(hour); 
        arcmin = deg / 60.0;
        arcsec = deg / 3600.0;

        Re = 6378137.0;
        wie = 7.2921151467e-5;
        f = 1.0 / 298.257223563;
        e2 = 2.0 * f - f * f;
        ge = 9.780325333434361; // Equatorial gravity
        g0 = 9.7803267715;      // Standard gravity reference

        m = (wie * wie * Re) / ge;
        beta = 2.5 * m - f;
        beta1 = 0.125 * (5.0 * m * f - f * f);

        ug = 1e-6 * g0; 
        Hz = 1.0;
        ugpsHz = ug / sqrt(Hz);
        ugpg2 = ug / (g0 * g0);
        nm = Re * arcmin;
    }
};