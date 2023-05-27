#include <stdio.h>
#include <math.h>
#include <stdlib.h>

struct Result{
    double value1;
    double value2;
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
struct Result vinc(double latp, double latc, double longp, double longc) {
    double req = 6378137.0;             //Radius at equator
    double flat = 1 / 298.257223563;    //flattenig of earth
    double rpol = (1 - flat) * req;

    double u1, u2, lon, lam, tol, diff, sin_sigma, cos_sigma, sigma, sin_alpha, cos_sq_alpha, cos2sigma;
    double A, B, C, lam_pre, delta_sig, dis, azi1, usq;

    // convert to radians
    latp = M_PI * latp / 180.0;
    latc = M_PI * latc / 180.0;
    longp = M_PI * longp / 180.0;
    longc = M_PI * longc / 180.0;

    u1 = atan((1 - flat) * tan(latc));
    u2 = atan((1 - flat) * tan(latp));

    lon = longp - longc;
    lam = lon;
    tol = pow(10., -12.); // iteration tolerance
    diff = 1.;

    while (fabs(diff) > tol) {
        sin_sigma = sqrt(pow((cos(u2) * sin(lam)), 2.) + pow(cos(u1)*sin(u2) - sin(u1)*cos(u2)*cos(lam), 2.));
        cos_sigma = sin(u1) * sin(u2) + cos(u1) * cos(u2) * cos(lam);
        sigma = atan(sin_sigma / cos_sigma);
        if (sigma <= 0) sigma = M_PI + sigma;
        sin_alpha = (cos(u1) * cos(u2) * sin(lam)) / sin_sigma;
        cos_sq_alpha = 1 - pow(sin_alpha, 2.);
        cos2sigma = cos_sigma - ((2 * sin(u1) * sin(u2)) / cos_sq_alpha);
        C = (flat / 16) * cos_sq_alpha * (4 + flat * (4 - 3 * cos_sq_alpha));
        lam_pre = lam;
        lam = lon + (1 - C) * flat * sin_alpha * (sigma + C * sin_sigma * (cos2sigma + C * cos_sigma * (2 * pow(cos2sigma, 2.) - 1)));
        diff = fabs(lam_pre - lam);
    }

    sin_sigma = sqrt(pow((cos(u2) * sin(lam)), 2.) + pow(cos(u1)*sin(u2) - sin(u1)*cos(u2)*cos(lam), 2.));
    cos_sigma = sin(u1) * sin(u2) + cos(u1) * cos(u2) * cos(lam);
    sigma = atan(sin_sigma / cos_sigma);
    if (sigma <= 0) sigma = M_PI + sigma;
    sin_alpha = (cos(u1) * cos(u2) * sin(lam)) / sin_sigma;
    cos_sq_alpha = 1 - pow(sin_alpha, 2.);
    cos2sigma = cos_sigma - ((2 * sin(u1) * sin(u2)) / cos_sq_alpha);

    usq = cos_sq_alpha * ((pow(req, 2.) - pow(rpol, 2.)) / pow(rpol ,2.));
    A = 1 + (usq / 16384) * (4096 + usq * (-768 + usq * (320 - 175 * usq)));
    B = (usq / 1024) * (256 + usq * (-128 + usq * (74 - 47 * usq)));
    delta_sig = B * sin_sigma * (cos2sigma + 0.25 * B * (cos_sigma * (-1 + 2 * pow(cos2sigma, 2.)) -
                                                         (1 / 6) * B * cos2sigma * (-3 + 4 * pow(sin_sigma, 2.)) *
                                                         (-3 + 4 * pow(cos2sigma, 2.))));
    dis = rpol * A * (sigma - delta_sig);
    azi1 = atan2((cos(u2) * sin(lam)), (cos(u1) * sin(u2) - sin(u1) * cos(u2) * cos(lam)));

    struct Result res = {dis, azi1};
    return(res);
}

/**
 * Transform longitude, latitude coordinates of two points to x,y coordinates
 * @param latp
 * @param latc
 * @param longp
 * @param longc
 * @return
 */
struct Result trans(double latp, double latc, double longp, double longc){

    double rav, theta, dis, azi, xy, x, y;
    struct Result dis_azi;

    rav = 6371000.0;  //average radius
    dis_azi = vinc(latp,latc,longp,longc);

    dis = dis_azi.value1;
    azi = dis_azi.value2;

    theta = dis/rav; //finding theta angle
    xy = sin(theta)*rav; //length in xy plane

    y = xy*cos(azi); //lat for chunk
    x = xy*sin(azi); //long for chunk
    struct Result res = {y, x};
    return(res);
}

int main() {
    // Testing with one value
    struct Result res = trans(40.99698, 46.0, 9.20127, 10.0);
    printf("%f - %f\n", res.value1, res.value2);
    return 0;
}
