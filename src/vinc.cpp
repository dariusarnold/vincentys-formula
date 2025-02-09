#include "vinc.hpp"

#include <tuple>
#include <cmath>

#ifdef VINC_PYBIND11
#include <pybind11/pybind11.h>
namespace py = pybind11;
#endif

std::tuple<double, double> vinc(double latp, double latc, double longp, double longc) {
    using namespace std;
    constexpr double req = 6378137.0;             //Radius at equator
    constexpr double flat = 1 / 298.257223563;    //flattening of earth
    constexpr double rpol = (1 - flat) * req;

    double sin_sigma, cos_sigma, sigma, sin_alpha, cos_sq_alpha, cos2sigma;
    double C, lam_pre;

    // convert to radians
    latp = M_PI * latp / 180.0;
    latc = M_PI * latc / 180.0;
    longp = M_PI * longp / 180.0;
    longc = M_PI * longc / 180.0;

    const double u1 = atan((1 - flat) * tan(latc));
    const double u2 = atan((1 - flat) * tan(latp));

    double lon = longp - longc;
    double lam = lon;
    constexpr double tol = 10.e-12; // iteration tolerance
    double diff = 1.;

    while (abs(diff) > tol) {
        sin_sigma = sqrt(pow((cos(u2) * sin(lam)), 2.) + pow(cos(u1)*sin(u2) - sin(u1)*cos(u2)*cos(lam), 2.));
        if (sin_sigma == 0.) {
            // Coincident points, prevent division by zero resulting in NaN.
            return {0., 0.};
        }
        cos_sigma = sin(u1) * sin(u2) + cos(u1) * cos(u2) * cos(lam);
        sigma = atan(sin_sigma / cos_sigma);
        if (sigma <= 0) sigma = M_PI + sigma;
        sin_alpha = (cos(u1) * cos(u2) * sin(lam)) / sin_sigma;
        cos_sq_alpha = 1 - pow(sin_alpha, 2.);
        if (cos_sq_alpha == 0.) {
            cos2sigma = 0.;
        } else {
            cos2sigma = cos_sigma - ((2 * sin(u1) * sin(u2)) / cos_sq_alpha);
        }
        C = (flat / 16) * cos_sq_alpha * (4 + flat * (4 - 3 * cos_sq_alpha));
        lam_pre = lam;
        lam = lon + (1 - C) * flat * sin_alpha * (sigma + C * sin_sigma * (cos2sigma + C * cos_sigma * (2 * pow(cos2sigma, 2.) - 1)));
        diff = abs(lam_pre - lam);
    }

    const double usq = cos_sq_alpha * ((pow(req, 2.) - pow(rpol, 2.)) / pow(rpol ,2.));
    const double A = 1 + (usq / 16384) * (4096 + usq * (-768 + usq * (320 - 175 * usq)));
    const double B = (usq / 1024) * (256 + usq * (-128 + usq * (74 - 47 * usq)));
    const double delta_sig = B * sin_sigma * (cos2sigma + 0.25 * B * (cos_sigma * (-1 + 2 * pow(cos2sigma, 2.)) -
                                                         (1. / 6) * B * cos2sigma * (-3 + 4 * pow(sin_sigma, 2.)) *
                                                         (-3 + 4 * pow(cos2sigma, 2.))));
    const double dis = rpol * A * (sigma - delta_sig);
    const double azi1 = atan2((cos(u2) * sin(lam)), (cos(u1) * sin(u2) - sin(u1) * cos(u2) * cos(lam)));

    return std::make_tuple(dis, azi1);
}

#ifdef VINC_PYBIND11
PYBIND11_MODULE(vinc_cpp, m){
    m.doc() = "Vincentys formula in C++";

    m.def("vinc", &vinc, "Vincentys formula");
}
#endif
