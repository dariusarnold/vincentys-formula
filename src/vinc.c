#include "vinc.h"

#include <stdio.h>
#include <math.h>


struct Geom {
    double sin_sigma;
    double cos_sigma;
    double sigma;
    double sin_alpha;
    double cos_sq_alpha;
    double cos2sigma;
};


void calcgeo(struct Geom *g, double lam, double u1, double u2) {
    g->sin_sigma = sqrt(pow((cos(u2) * sin(lam)), 2.) + pow(cos(u1) * sin(u2) - sin(u1) * cos(u2) * cos(lam), 2.));
    g->cos_sigma = sin(u1) * sin(u2) + cos(u1) * cos(u2) * cos(lam);
    g->sigma = atan(g->sin_sigma / g->cos_sigma);
    if (g->sigma <= 0) {
        g->sigma += M_PI;
    }
    g->sin_alpha = (cos(u1) * cos(u2) * sin(lam)) / g->sin_sigma;
    g->cos_sq_alpha = 1 - pow(g->sin_alpha, 2.);
    g->cos2sigma = g->cos_sigma - ((2 * sin(u1) * sin(u2)) / g->cos_sq_alpha);
}


struct ResultVinc vinc(double latp, double latc, double longp, double longc) {
    struct Geom gvar;
    double req = 6378137.0;             // Radius at equator
    double flat = 1 / 298.257223563;    // Flattening of earth
    double rpol = (1 - flat) * req;

    double u1, u2, lon, lam, tol, diff;
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
    tol = pow(10., -12.); // Iteration tolerance, 10E-12 corresponds to approx. 0.06 mm
    diff = 1.;

    while (fabs(diff) > tol) {
        calcgeo(&gvar, lam, u1, u2);
        C = (flat / 16) * gvar.cos_sq_alpha * (4 + flat * (4 - 3 * gvar.cos_sq_alpha));
        lam_pre = lam;
        lam = lon + (1 - C) * flat * gvar.sin_alpha * (gvar.sigma + C * gvar.sin_sigma * (gvar.cos2sigma + C * gvar.cos_sigma * (2 * pow(gvar.cos2sigma, 2.) - 1)));
        diff = fabs(lam_pre - lam);
    }

    calcgeo(&gvar, lam, u1, u2);

    usq = gvar.cos_sq_alpha * ((pow(req, 2.) - pow(rpol, 2.)) / pow(rpol ,2.));
    A = 1 + (usq / 16384) * (4096 + usq * (-768 + usq * (320 - 175 * usq)));
    B = (usq / 1024) * (256 + usq * (-128 + usq * (74 - 47 * usq)));
    delta_sig = B * gvar.sin_sigma * (gvar.cos2sigma + 0.25 * B * (gvar.cos_sigma * (-1 + 2 * pow(gvar.cos2sigma, 2.)) -
                                                                   (1 / 6) * B * gvar.cos2sigma * (-3 + 4 * pow(gvar.sin_sigma, 2.)) *
                                                                   (-3 + 4 * pow(gvar.cos2sigma, 2.))));
    dis = rpol * A * (gvar.sigma - delta_sig);
    azi1 = atan2((cos(u2) * sin(lam)), (cos(u1) * sin(u2) - sin(u1) * cos(u2) * cos(lam)));

    struct ResultVinc r = {dis, azi1};
    return r;
}


struct ResultTrans trans(double latp, double latc, double longp, double longc){

    double rav, theta, xy, x, y;
    struct ResultVinc dis_azi;

    rav = 6371000.0;  // average radius
    dis_azi = vinc(latp,latc,longp,longc);

    theta = dis_azi.distance / rav; // finding theta angle
    xy = sin(theta) * rav; // length in xy plane
    y = xy * cos(dis_azi.distance); // lat for chunk
    x = xy * sin(dis_azi.azimuth); // long for chunk
    struct ResultTrans r = {x, y};
    return r;
}
