# topocent.m
import numpy as np


# togeod.m
def togeod(a=None, finv=None, X=None, Y=None, Z=None, *args, **kwargs):
    # TOGEOD   Subroutine to calculate geodetic coordinates latitude, longitude,
    #         height given Cartesian coordinates X,Y,Z, and reference ellipsoid
    #         values semi-major axis (a) and the inverse of flattening (finv).

    # [dphi, dlambda, h] = togeod(a, finv, X, Y, Z);

    #  The units of linear parameters X,Y,Z,a must all agree (m,km,mi,ft,..etc)
    #  The output units of angular quantities will be in decimal degrees
    #  (15.5 degrees not 15 deg 30 min). The output units of h will be the
    #  same as the units of X,Y,Z,a.

    #   Inputs:
    #       a           - semi-major axis of the reference ellipsoid
    #       finv        - inverse of flattening of the reference ellipsoid
    #       X,Y,Z       - Cartesian coordinates

    #   Outputs:
    #       dphi        - latitude
    #       dlambda     - longitude
    #       h           - height above reference ellipsoid

    #  Copyright (C) 1987 C. Goad, Columbus, Ohio
    #  Reprinted with permission of author, 1996
    #  Fortran code translated into MATLAB
    #  Kai Borre 03-30-96

    # CVS record:
    # $Id: togeod.m,v 1.1.1.1.2.4 2006/08/22 13:45:59 dpl Exp $
    # ==========================================================================

    h = 0.0
    # togeod.m:32
    tolsq = 1e-10
    # togeod.m:33
    maxit = 10
    # togeod.m:34
    # compute radians-to-degree factor
    rtd = 180 / np.pi
    # togeod.m:37
    # compute square of eccentricity
    if finv < 1e-20:
        esq = 0.0
    # togeod.m:41
    else:
        esq = (2 - 1 / finv) / finv
    # togeod.m:43

    oneesq = 1 - esq
    # togeod.m:46
    # first guess
    # P is distance from spin axis
    P = np.sqrt(X ** 2 + Y ** 2)
    # togeod.m:50
    # direct calculation of longitude

    if P > 1e-20:
        dlambda = np.arctan2(Y, X) * rtd
    # togeod.m:54
    else:
        dlambda = 0.0
    # togeod.m:56

    if dlambda < 0:
        dlambda = dlambda + 360
    # togeod.m:60

    # r is distance from origin (0,0,0)
    r = np.sqrt(P ** 2 + Z ** 2)
    # togeod.m:64
    if r > 1e-20:
        sinphi = Z / r
    # togeod.m:67
    else:
        sinphi = 0.0
    # togeod.m:69

    dphi = np.arcsin(sinphi)
    # togeod.m:72
    # initial value of height  =  distance from origin minus
    # approximate distance from origin to surface of ellipsoid
    if r < 1e-20:
        h = 0.0
        # togeod.m:77
        return dphi, dlambda, h

    h = r - a * (1 - sinphi * sinphi / finv)
    # togeod.m:81
    # iterate
    for i in range(maxit):
        sinphi = np.sin(dphi)
        # togeod.m:85
        cosphi = np.cos(dphi)
        # togeod.m:86
        N_phi = a / np.sqrt(1 - esq * sinphi * sinphi)
        # togeod.m:89
        dP = P - (N_phi + h) * cosphi
        # togeod.m:92
        dZ = Z - (N_phi * oneesq + h) * sinphi
        # togeod.m:93
        h = h + sinphi * dZ + cosphi * dP
        # togeod.m:96
        dphi = dphi + (cosphi * dZ - sinphi * dP) / (N_phi + h)
        # togeod.m:97
        if (dP * dP + dZ * dZ) < tolsq:
            break
        # Not Converged--Warn user
        if i == maxit - 1:
            print ' Problem in TOGEOD, did not converge in %2.0f iterations' % i

    dphi *= rtd
    return dphi, dlambda, h


# togeod.m:111
######## end togeod.m  ######################


def topocent(X=None, dx=None, *args, **kwargs):
    # TOPOCENT  Transformation of vector dx into topocentric coordinate
    #          system with origin at X.
    #          Both parameters are 3 by 1 vectors.

    # [Az, El, D] = topocent(X, dx);

    #   Inputs:
    #       X           - vector origin corrdinates (in ECEF system [X; Y; Z;])
    #       dx          - vector ([dX; dY; dZ;]).

    #   Outputs:
    #       D           - vector length. Units like units of the input
    #       Az          - azimuth from north positive clockwise, degrees
    #       El          - elevation angle, degrees

    # Kai Borre 11-24-96
    # Copyright (c) by Kai Borre

    # CVS record:
    # $Id: topocent.m,v 1.1.1.1.2.4 2006/08/22 13:45:59 dpl Exp $
    # ==========================================================================

    dtr = np.pi / 180
    # topocent.m:24
    phi, lambda_, h = togeod(6378137, 298.257223563, X[0], X[1], X[2])
    # topocent.m:26
    cl = np.cos(lambda_ * dtr)
    # topocent.m:28
    sl = np.sin(lambda_ * dtr)
    # topocent.m:29
    cb = np.cos(phi * dtr)
    # topocent.m:30
    sb = np.sin(phi * dtr)
    # topocent.m:31
    F = np.array([[- sl, -sb * cl, cb * cl], [cl, -sb * sl, cb * sl], [0.0, cb, sb]])
    # topocent.m:33
    local_vector = F.T.dot(dx)
    # topocent.m:37
    E = local_vector[0]
    # topocent.m:38
    N = local_vector[1]
    # topocent.m:39
    U = local_vector[2]
    # topocent.m:40
    hor_dis = np.sqrt(E ** 2 + N ** 2)
    # topocent.m:42
    if hor_dis < 1e-20:
        Az = 0.0
        # topocent.m:45
        El = 90.0
    # topocent.m:46
    else:
        Az = np.arctan2(E, N) / dtr
        # topocent.m:48
        El = np.arctan2(U, hor_dis) / dtr
    # topocent.m:49

    if Az < 0:
        Az = Az + 360
    # topocent.m:53

    D = np.sqrt(dx[0] ** 2 + dx[1] ** 2 + dx[2] ** 2)
    return Az, El, D
# topocent.m:56
######### end topocent.m #########
