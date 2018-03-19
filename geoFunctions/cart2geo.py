# cart2geo.m
import numpy as np


def cart2geo(X=None, Y=None, Z=None, i=None, *args, **kwargs):
    # CART2GEO Conversion of Cartesian coordinates (X,Y,Z) to geographical
    # coordinates (phi, lambda, h) on a selected reference ellipsoid.

    # [phi, lambda, h] = cart2geo(X, Y, Z, i);

    #   Choices i of Reference Ellipsoid for Geographical Coordinates
    #          	  1. International Ellipsoid 1924
    #	          2. International Ellipsoid 1967
    #	          3. World Geodetic System 1972
    #	          4. Geodetic Reference System 1980
    #	          5. World Geodetic System 1984

    # Kai Borre 10-13-98
    # Copyright (c) by Kai Borre
    # Revision: 1.0   Date: 1998/10/23

    # CVS record:
    # $Id: cart2geo.m,v 1.1.2.3 2007/01/29 15:22:49 dpl Exp $
    # ==========================================================================

    a = np.array([6378388.0, 6378160.0, 6378135.0, 6378137.0, 6378137.0])
    # cart2geo.m:22
    f = np.array([1 / 297, 1 / 298.247, 1 / 298.26, 1 / 298.257222101, 1 / 298.257223563])
    # cart2geo.m:23
    lambda_ = np.arctan2(Y, X)
    # cart2geo.m:25
    ex2 = (2 - f[i]) * f[i] / ((1 - f[i]) ** 2)
    # cart2geo.m:26
    c = a[i] * np.sqrt(1 + ex2)
    # cart2geo.m:27
    phi = np.arctan(Z / (np.sqrt(X ** 2 + Y ** 2) * (1 - (2 - f[i])) * f[i]))
    # cart2geo.m:28
    h = 0.1
    # cart2geo.m:30
    oldh = 0
    # cart2geo.m:30
    iterations = 0
    # cart2geo.m:31
    while abs(h - oldh) > 1e-12:

        oldh = h
        # cart2geo.m:33
        N = c / np.sqrt(1 + ex2 * np.cos(phi) ** 2)
        # cart2geo.m:34
        phi = np.arctan(Z / (np.sqrt(X ** 2 + Y ** 2) * (1 - (2 - f[i]) * f[i] * N / (N + h))))
        # cart2geo.m:35
        h = np.sqrt(X ** 2 + Y ** 2) / np.cos(phi) - N
        # cart2geo.m:36
        iterations += 1
        # cart2geo.m:38
        if iterations > 100:
            print 'Failed to approximate h with desired precision. h-oldh: %e.' % (h - oldh)
            break

    phi *= (180 / np.pi)
    # cart2geo.m:45
    # b = zeros(1,3);
    # b(1,1) = fix(phi);
    # b(2,1) = fix(rem(phi,b(1,1))*60);
    # b(3,1) = (phi-b(1,1)-b(1,2)/60)*3600;

    lambda_ *= (180 / np.pi)
    # cart2geo.m:51
    # l = zeros(1,3);
    # l(1,1) = fix(lambda);
    # l(2,1) = fix(rem(lambda,l(1,1))*60);
    # l(3,1) = (lambda-l(1,1)-l(1,2)/60)*3600;

    # fprintf('\n     phi =#3.0f #3.0f #8.5f',b(1),b(2),b(3))
    # fprintf('\n  lambda =#3.0f #3.0f #8.5f',l(1),l(2),l(3))
    # fprintf('\n       h =#14.3f\n',h)
    return phi, lambda_, h
############## end cart2geo.m ###################
