#pragma once


struct ResultVinc {
    double distance;
    double azimuth;
};


struct ResultTrans {
    double x;
    double y;
};


/**
 * Vincenty's formulae for distance, solving the inverse problem.
 * Calculates distance and azimuth between two points on the surface of a spheroid.
 * For more information see: https://en.wikipedia.org/wiki/Vincenty%27s_formulae
 * @param latp Latitude
 * @param latc Latitude
 * @param longp Longitude
 * @param longc Longitude
 * @return
 */
struct ResultVinc vinc(double latp, double latc, double longp, double longc);


/**
 * Transform longitude, latitude coordinates of two points to x,y coordinates
 * @param latp
 * @param latc
 * @param longp
 * @param longc
 * @return
 */
struct ResultTrans trans(double latp, double latc, double longp, double longc);
