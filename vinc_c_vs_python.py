#!/usr/bin/env python3

import argparse
import ctypes
import pathlib
import sys
import types
from ctypes import c_double
import numpy as np
import timeit
from math import sin, cos, atan, atan2, tan, sqrt


# vinc function uses a struct containing azimut and distance as members, recreate this
class Result(ctypes.Structure):
    _fields_ = [("value1", c_double), ("value2", c_double)]


# vinc numpy math functions
def vinc_numpy(latp: float, latc: float, longp: float, longc: float) -> tuple[float, float]:

    req = 6378137.0 #Radius at equator
    flat = 1/298.257223563 #flattenig of earth
    rpol = (1-flat)*req

    #convert to radians
    latp = np.pi*latp/180.0
    latc = np.pi*latc/180.0
    longp = np.pi*longp/180.0
    longc = np.pi*longc/180.0

    u1 = np.arctan((1-flat)*np.tan(latc))
    u2 = np.arctan((1-flat)*np.tan(latp))

    lon = longp-longc

    lam = lon

    tol = 10**(-12) # iteration tolerance
    max_iterations = 1_000_000

    for _ in range(max_iterations):
        sin_sigma = np.sqrt((np.cos(u2)*np.sin(lam))**2 + (np.cos(u1)*np.sin(u2) - np.sin(u1)*np.cos(u2)*np.cos(lam))**2)
        cos_sigma = np.sin(u1)*np.sin(u2) + np.cos(u1)*np.cos(u2)*np.cos(lam)
        sigma = np.arctan(sin_sigma/cos_sigma)
        sin_alpha = (np.cos(u1)*np.cos(u2)*np.sin(lam))/sin_sigma
        cos_sq_alpha = 1- sin_alpha**2
        cos2sigma = cos_sigma -((2*np.sin(u1)*np.sin(u2))/cos_sq_alpha)
        C = (flat/16)*cos_sq_alpha*(4 + flat*(4 - 3*cos_sq_alpha))
        lam_pre = lam
        lam = lon + (1-C)*flat*sin_alpha*(sigma+C*sin_sigma*(cos2sigma + C*cos_sigma*(2*cos2sigma**2 - 1)))

        diff = abs(lam_pre - lam)
        if abs(diff) < tol:
            break
    else:
        raise RuntimeError(f"Tolerance value chosen to large ({tol}).")

    usq = cos_sq_alpha*((req**2-rpol**2)/rpol**2)
    A = 1 + (usq/16384)*(4096+usq*(-768+usq*(320-175*usq)))
    B = (usq/1024)*(256+usq*(-128+usq*(74-47*usq)))
    delta_sig = B*sin_sigma*(cos2sigma+0.25*B*(cos_sigma*(-1+2*cos2sigma**2)-(1/6)*B*cos2sigma*(-3+4*sin_sigma**2)*(-3+4*cos2sigma**2)))

    dis = rpol*A*(sigma-delta_sig)

    azi1 = np.arctan2((np.cos(u2)*np.sin(lam)), (np.cos(u1)*np.sin(u2) - np.sin(u1)*np.cos(u2)*np.cos(lam)))
    #azi2 = np.arctan2((np.cos(u1)*np.sin(lam))/(-np.sin(u1)*np.cos(u2) + np.cos(u1)*np.sin(u2)*np.cos(lam)))

    return dis,azi1


# python implementation of vinc
def vinc_pure_python(latp: float, latc: float, longp: float, longc: float) -> tuple[float, float]:

    req = 6378137.0 #Radius at equator
    flat = 1/298.257223563 #flattenig of earth
    rpol = (1-flat)*req

    #convert to radians
    latp = np.pi*latp/180.0
    latc = np.pi*latc/180.0
    longp = np.pi*longp/180.0
    longc = np.pi*longc/180.0

    u1 = atan((1-flat)*tan(latc))
    u2 = atan((1-flat)*tan(latp))

    lon = longp-longc

    lam = lon

    tol = 10**(-12) # iteration tolerance
    max_iterations = 1_000_000

    for _ in range(max_iterations):
        sin_sigma = sqrt((cos(u2)*sin(lam))**2 + (cos(u1)*sin(u2) - sin(u1)*cos(u2)*cos(lam))**2)
        cos_sigma = sin(u1)*sin(u2) + cos(u1)*cos(u2)*cos(lam)
        sigma = atan(sin_sigma/cos_sigma)
        sin_alpha = (cos(u1)*cos(u2)*sin(lam))/sin_sigma
        cos_sq_alpha = 1- sin_alpha**2
        cos2sigma = cos_sigma -((2*sin(u1)*sin(u2))/cos_sq_alpha)
        C = (flat/16)*cos_sq_alpha*(4 + flat*(4 - 3*cos_sq_alpha))
        lam_pre = lam
        lam = lon + (1-C)*flat*sin_alpha*(sigma+C*sin_sigma*(cos2sigma + C*cos_sigma*(2*cos2sigma**2 - 1)))

        diff = abs(lam_pre - lam)
        if abs(diff) < tol:
            break
    else:
        raise RuntimeError(f"Tolerance value chosen to large ({tol}).")

    usq = cos_sq_alpha*((req**2-rpol**2)/rpol**2)
    A = 1 + (usq/16384)*(4096+usq*(-768+usq*(320-175*usq)))
    B = (usq/1024)*(256+usq*(-128+usq*(74-47*usq)))
    delta_sig = B*sin_sigma*(cos2sigma+0.25*B*(cos_sigma*(-1+2*cos2sigma**2)-(1/6)*B*cos2sigma*(-3+4*sin_sigma**2)*(-3+4*cos2sigma**2)))

    dis = rpol*A*(sigma-delta_sig)

    azi1 = atan2((cos(u2)*sin(lam)), (cos(u1)*sin(u2) - sin(u1)*cos(u2)*cos(lam)))
    #azi2 = np.arctan2((np.cos(u1)*np.sin(lam))/(-np.sin(u1)*np.cos(u2) + np.cos(u1)*np.sin(u2)*np.cos(lam)))

    return dis,azi1


def trans(latp: float, latc: float, longp:float, longc:float) -> tuple[float, float]:

    rav = 6371000.0 #average radius

    dis, azi = vinc_numpy(latp,latc,longp,longc)

    theta = dis/rav #finding theta angle
    xy = np.sin(theta)*rav #length in xy plane

    y = xy*np.cos(azi) #lat for chunk
    x = xy*np.sin(azi) #long for chunk

    return y, x


def test_vinc(vincer: ctypes.CDLL, vinc_cpp: types.ModuleType) -> None:
    number_of_calls = 1000

    print("\nTesting function vinc")

    # setup for python vinc function
    setup_str_python = "latp, latc, longp, longc = 40.99698, 46.0, 9.20127, 10.0; from __main__ import vinc_pure_python"
    # python code which is called
    stmt_python = 'vinc_pure_python(latp, latc, longp, longc)'

    # setup for numpy vinc function
    setup_str_numpy = "latp, latc, longp, longc = 40.99698, 46.0, 9.20127, 10.0; from __main__ import vinc_numpy"
    # python code which is called
    stmt_numpy = 'vinc_numpy(latp, latc, longp, longc)'

    # setup for c vinc function
    setup_str_c = "latp, latc, longp, longc = 40.99698, 46.0, 9.20127, 10.0"
    # python code which calls c vinc function
    stmt_c = 'vincer.vinc(c_double(latp), c_double(latc), c_double(longp), c_double(longc))'

    # setup for the cpp vinc function
    setup_str_cpp = "latp, latc, longp, longc = 40.99698, 46.0, 9.20127, 10.0"
    # python code to call cpp function
    stmt_cpp = "vinc_cpp.vinc(latp, latc, longp, longc)"

    times_python = timeit.repeat(stmt=stmt_python, setup=setup_str_python, repeat=10, number=number_of_calls)
    times_numpy = timeit.repeat(stmt=stmt_numpy, setup=setup_str_numpy, repeat=10, number=number_of_calls)
    times_c = timeit.repeat(stmt=stmt_c, setup=setup_str_c, repeat=10, number=number_of_calls,
                            globals={"vincer": vincer, "c_double": c_double})
    times_cpp = timeit.repeat(stmt=stmt_cpp, setup=setup_str_cpp, repeat=10, number=number_of_calls,
                              globals={"vinc_cpp": vinc_cpp})

    # calculate time per run
    times_python = [res / number_of_calls for res in times_python]
    times_numpy = [res / number_of_calls for res in times_numpy]
    times_c = [res / number_of_calls for res in times_c]
    times_cpp = [res / number_of_calls for res in times_cpp]

    # Calculate results to compare numerical accuracy
    latp, latc, longp, longc = 40.99698, 46.0, 9.20127, 10.0
    result_python = vinc_pure_python(latp, latc, longp, longc)
    result_numpy = vinc_numpy(latp, latc, longp, longc)
    result_c = vincer.vinc(c_double(latp), c_double(latc), c_double(longp), c_double(longc))
    result_cpp = vinc_cpp.vinc(latp, latc, longp, longc)

    # print time and numerical results
    print_results(times_python, times_numpy, times_c, times_cpp, result_python, result_numpy, result_c, result_cpp, number_of_calls)

def test_trans(vincer: ctypes.CDLL) -> None:
    raise RuntimeError("Currently not supported")
    number_of_calls = 1000

    print("Testing function trans")

    # setup for python trans function
    setup_str_python = "latp, latc, longp, longc = 40.99698, 46.0, 9.20127, 10.0; from __main__ import trans"
    # python code which is called
    stmt_python = 'trans(latp, latc, longp, longc)'

    # setup for c vinc function
    setup_str_c = "latp, latc, longp, longc = 40.99698, 46.0, 9.20127, 10.0; from __main__ import vincer; from ctypes import c_double"
    # python code which calls c vinc function
    stmt_c = 'vincer.trans(c_double(latp), c_double(latc), c_double(longp), c_double(longc))'

    times_python = timeit.repeat(stmt=stmt_python, setup=setup_str_python, repeat=3, number=number_of_calls)
    times_c = timeit.repeat(stmt=stmt_c, setup=setup_str_c, repeat=3, number=number_of_calls)

    # calculate time per run
    times_python = [res / number_of_calls for res in times_python]
    times_c = [res / number_of_calls for res in times_c]

    # Calculate results to compare numerical accuracy
    latp, latc, longp, longc = 40.99698, 46.0, 9.20127, 10.0
    result_python = trans(latp, latc, longp, longc)
    result_c = vincer.trans(c_double(latp), c_double(latc), c_double(longp), c_double(longc))

    # print time and numerical results
    print_results(times_python, times_c, result_python, result_c, number_of_calls)


def load_vinc_so(so_path: pathlib.Path) -> ctypes.CDLL:
    # load the library containing the vinc function
    so_path = so_path.resolve()
    print(f"Trying to load {so_path}")
    vincer = ctypes.CDLL(str(so_path))
    # set the expected return type of the function vinc to result
    vincer.vinc.restype = Result
    vincer.trans.restype = Result
    return vincer


def print_results(times_python: list[float], times_numpy: list[float], times_c: list[float], times_cpp: list[float],
                  result_python: tuple[float, float], result_numpy: tuple[float, float], result_c: tuple[float, float],
                  result_cpp: Result,
                  number_of_calls: int) -> None:
    print("NumPy best average out of {} runs:  {:.1E} secs".format(number_of_calls, min(times_numpy)))
    print("Python best average out of {} runs: {:.1E} secs".format(number_of_calls, min(times_python)))
    print("C best average out of {} runs:      {:.1E} secs".format(number_of_calls, min(times_c)))
    print("C++ best average out of {} runs:    {:.1E} secs".format(number_of_calls, min(times_cpp)))
    print()
    print("Numpy results:  {} {}".format(*result_numpy))
    print("Python results: {} {}".format(*result_python))
    print("C results:      {} {}".format(result_c.value1, result_c.value2))
    print(f"C++ results:    {result_cpp[0]} {result_cpp[1]}")


def main():
    parser = argparse.ArgumentParser(
        prog="VincentyCompare",
        description="Vincenty Comparison Python vs C"
    )

    parser.add_argument("VINC_SO", type=pathlib.Path,
                        help="Path to vinc shared library file to use, typically libvinc.so.")
    parser.add_argument("VINC_CPP", type=pathlib.Path,
                        help="Path to C++ Python extension module, for example vinc_cpp.cpython-311-x86_64-linux-gnu.so.")
    args = parser.parse_args()
    vincer = load_vinc_so(args.VINC_SO)
    vinc_cpp_dir: pathlib.Path = args.VINC_CPP.absolute().parent
    print(f"Adding {vinc_cpp_dir} to sys.path")
    sys.path.append(str(vinc_cpp_dir))
    import vinc_cpp

    test_vinc(vincer, vinc_cpp)
    #test_trans(vincer)


if __name__ == '__main__':
    main()
