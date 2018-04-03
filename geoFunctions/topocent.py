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

    # ==========================================================================

    h = 0.0

    tolsq = 1e-10

    maxit = 10

    # compute radians-to-degree factor
    rtd = 180 / np.pi

    # compute square of eccentricity
    if finv < 1e-20:
        esq = 0.0

    else:
        esq = (2 - 1 / finv) / finv

    oneesq = 1 - esq

    # first guess
    # P is distance from spin axis
    P = np.sqrt(X ** 2 + Y ** 2)

    # direct calculation of longitude

    if P > 1e-20:
        dlambda = np.arctan2(Y, X) * rtd

    else:
        dlambda = 0.0

    if dlambda < 0:
        dlambda = dlambda + 360

    # r is distance from origin (0,0,0)
    r = np.sqrt(P ** 2 + Z ** 2)

    if r > 1e-20:
        sinphi = Z / r

    else:
        sinphi = 0.0

    dphi = np.arcsin(sinphi)

    # initial value of height  =  distance from origin minus
    # approximate distance from origin to surface of ellipsoid
    if r < 1e-20:
        h = 0.0

        return dphi, dlambda, h

    h = r - a * (1 - sinphi * sinphi / finv)

    # iterate
    for i in range(maxit):
        sinphi = np.sin(dphi)

        cosphi = np.cos(dphi)

        N_phi = a / np.sqrt(1 - esq * sinphi * sinphi)

        dP = P - (N_phi + h) * cosphi

        dZ = Z - (N_phi * oneesq + h) * sinphi

        h = h + sinphi * dZ + cosphi * dP

        dphi = dphi + (cosphi * dZ - sinphi * dP) / (N_phi + h)

        if (dP * dP + dZ * dZ) < tolsq:
            break
        # Not Converged--Warn user
        if i == maxit - 1:
            print ' Problem in TOGEOD, did not converge in %2.0f iterations' % i

    dphi *= rtd
    return dphi, dlambda, h


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

    # ==========================================================================

    dtr = np.pi / 180

    phi, lambda_, h = togeod(6378137, 298.257223563, X[0], X[1], X[2])

    cl = np.cos(lambda_ * dtr)

    sl = np.sin(lambda_ * dtr)

    cb = np.cos(phi * dtr)

    sb = np.sin(phi * dtr)

    F = np.array([[- sl, -sb * cl, cb * cl], [cl, -sb * sl, cb * sl], [0.0, cb, sb]])

    local_vector = F.T.dot(dx)

    E = local_vector[0]

    N = local_vector[1]

    U = local_vector[2]

    hor_dis = np.sqrt(E ** 2 + N ** 2)

    if hor_dis < 1e-20:
        Az = 0.0

        El = 90.0

    else:
        Az = np.arctan2(E, N) / dtr

        El = np.arctan2(U, hor_dis) / dtr

    if Az < 0:
        Az = Az + 360

    D = np.sqrt(dx[0] ** 2 + dx[1] ** 2 + dx[2] ** 2)
    return Az, El, D

######### end topocent.m #########
