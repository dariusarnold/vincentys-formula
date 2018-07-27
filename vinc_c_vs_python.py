#!/usr/bin/env python

from ctypes import *
import numpy as np
import timeit
from math import sin, cos, atan, atan2, tan, sqrt


# load the library containing the vinc function
vincer = CDLL("./vinc.so")

# vinc function uses a struct containing azimut and distance as members, recreate this
class Result(Structure):
    _fields_ = [("value1", c_double), ("value2", c_double)]
# set the expected return type of the function vinc to result
vincer.vinc.restype = Result
vincer.trans.restype = Result

# vinc numpy math functions
def vinc_numpy(latp, latc, longp, longc):

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

    tol = 10**(-12) # iteration tool
    diff = 1

    while abs(diff) > tol:

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

    usq = cos_sq_alpha*((req**2-rpol**2)/rpol**2)
    A = 1 + (usq/16384)*(4096+usq*(-768+usq*(320-175*usq)))
    B = (usq/1024)*(256+usq*(-128+usq*(74-47*usq)))
    delta_sig = B*sin_sigma*(cos2sigma+0.25*B*(cos_sigma*(-1+2*cos2sigma**2)-(1/6)*B*cos2sigma*(-3+4*sin_sigma**2)*(-3+4*cos2sigma**2)))

    dis = rpol*A*(sigma-delta_sig)

    azi1 = np.arctan2((np.cos(u2)*np.sin(lam)), (np.cos(u1)*np.sin(u2) - np.sin(u1)*np.cos(u2)*np.cos(lam)))
    #azi2 = np.arctan2((np.cos(u1)*np.sin(lam))/(-np.sin(u1)*np.cos(u2) + np.cos(u1)*np.sin(u2)*np.cos(lam)))

    return dis,azi1


# python implementation of vinc
def vinc_pure_python(latp, latc, longp, longc):

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

    tol = 10**(-12) # iteration tool
    diff = 1

    while abs(diff) > tol:

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

    usq = cos_sq_alpha*((req**2-rpol**2)/rpol**2)
    A = 1 + (usq/16384)*(4096+usq*(-768+usq*(320-175*usq)))
    B = (usq/1024)*(256+usq*(-128+usq*(74-47*usq)))
    delta_sig = B*sin_sigma*(cos2sigma+0.25*B*(cos_sigma*(-1+2*cos2sigma**2)-(1/6)*B*cos2sigma*(-3+4*sin_sigma**2)*(-3+4*cos2sigma**2)))

    dis = rpol*A*(sigma-delta_sig)

    azi1 = atan2((cos(u2)*sin(lam)), (cos(u1)*sin(u2) - sin(u1)*cos(u2)*cos(lam)))
    #azi2 = np.arctan2((np.cos(u1)*np.sin(lam))/(-np.sin(u1)*np.cos(u2) + np.cos(u1)*np.sin(u2)*np.cos(lam)))

    return dis,azi1

def trans(latp,latc,longp,longc):

    rav = 6371000.0 #average radius

    dis, azi = vinc_numpy(latp,latc,longp,longc)

    theta = dis/rav #finding theta angle
    xy = np.sin(theta)*rav #length in xy plane

    y = xy*np.cos(azi) #lat for chunk
    x = xy*np.sin(azi) #long for chunk

    return y, x

def test_vinc_numpy_math():
    number_of_calls = 1000

    print("\nTesting function vinc with Numpy Math")

    # setup for python vinc function
    setup_str_python = "latp, latc, longp, longc = 40.99698, 46.0, 9.20127, 10.0; from __main__ import vinc_numpy"
    # python code which is called
    stmt_python = 'vinc_numpy(latp, latc, longp, longc)'

    # setup for c vinc function
    setup_str_c = "latp, latc, longp, longc = 40.99698, 46.0, 9.20127, 10.0; from __main__ import vincer; from ctypes import c_double"
    # python code which calls c vinc function
    stmt_c = 'vincer.vinc(c_double(latp), c_double(latc), c_double(longp), c_double(longc))'

    times_python = timeit.repeat(stmt=stmt_python, setup=setup_str_python, repeat=3, number=number_of_calls)
    times_c = timeit.repeat(stmt=stmt_c, setup=setup_str_c, repeat=3, number=number_of_calls)

    # calculate time per run
    times_python = [res / number_of_calls for res in times_python]
    times_c = [res / number_of_calls for res in times_c]

    # Calculate results to compare numerical accuracy
    latp, latc, longp, longc = 40.99698, 46.0, 9.20127, 10.0
    result_python = vinc_pure_python(latp, latc, longp, longc)
    result_c = vincer.vinc(c_double(latp), c_double(latc), c_double(longp), c_double(longc))

    # print time and numrical results
    print_results(times_python, times_c, result_python, result_c, number_of_calls)

def test_vinc_python_math():
    number_of_calls = 1000

    print("\nTesting function vinc with Python Math")

    # setup for python vinc function
    setup_str_python = "latp, latc, longp, longc = 40.99698, 46.0, 9.20127, 10.0; from __main__ import vinc_pure_python"
    # python code which is called
    stmt_python = 'vinc_pure_python(latp, latc, longp, longc)'

    # setup for c vinc function
    setup_str_c = "latp, latc, longp, longc = 40.99698, 46.0, 9.20127, 10.0; from __main__ import vincer; from ctypes import c_double"
    # python code which calls c vinc function
    stmt_c = 'vincer.vinc(c_double(latp), c_double(latc), c_double(longp), c_double(longc))'

    times_python = timeit.repeat(stmt=stmt_python, setup=setup_str_python, repeat=3, number=number_of_calls)
    times_c = timeit.repeat(stmt=stmt_c, setup=setup_str_c, repeat=3, number=number_of_calls)

    # calculate time per run
    times_python = [res / number_of_calls for res in times_python]
    times_c = [res / number_of_calls for res in times_c]

    # Calculate results to compare numerical accuracy
    latp, latc, longp, longc = 40.99698, 46.0, 9.20127, 10.0
    result_python = vinc_numpy(latp, latc, longp, longc)
    result_c = vincer.vinc(c_double(latp), c_double(latc), c_double(longp), c_double(longc))

    # print time and numrical results
    print_results(times_python, times_c, result_python, result_c, number_of_calls)

def test_trans():
    number_of_calls = 1000

    print("\nTesting function trans")

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

    # print time and numrical results
    print_results(times_python, times_c, result_python, result_c, number_of_calls)

def print_results(times_python, times_c, result_python, result_c, number_of_calls):
    print("Python best average out of {} runs: {:.12E} secs".format(number_of_calls, min(times_python)))
    print("C best average out of {} runs:      {:.12E} secs".format(number_of_calls, min(times_c)))
    print("C is {:4.2f} times faster".format(min(times_python) / min(times_c)))

    print("Python results: {} {}".format(*result_python))
    print("C results:      {} {}".format(result_c.value1, result_c.value2))
    print("Difference distance: {:.8E}%".format((result_python[0] / result_c.value1) * 100 - 100))
    print("Difference azimut:   {:.8E}%".format((result_python[1] / result_c.value2) * 100 - 100))

if __name__ == '__main__':
    test_vinc_numpy_math()
    test_vinc_python_math()
    test_trans()
