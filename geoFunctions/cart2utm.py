import numpy as np


# clsin.m
def clsin(ar=None, degree=None, argument=None, *args, **kwargs):
    # Clenshaw summation of sinus of argument.

    # result = clsin(ar, degree, argument);

    # Written by Kai Borre
    # December 20, 1995

    # See also WGS2UTM or CART2UTM

    # CVS record:
    # $Id: clsin.m,v 1.1.1.1.2.4 2006/08/22 13:45:59 dpl Exp $
    # ==========================================================================

    cos_arg = 2 * np.cos(argument)
    # clsin.m:15
    hr1 = 0
    # clsin.m:16
    hr = 0
    # clsin.m:17
    # TODO fix index range of t
    for t in range(degree, 0, -1):
        hr2 = hr1
        # clsin.m:20
        hr1 = hr
        # clsin.m:21
        hr = ar[t - 1] + cos_arg * hr1 - hr2
    # clsin.m:22

    result = hr * np.sin(argument)
    return result


# clsin.m:25
####################### end clsin.m  #####################


# clksin.m
def clksin(ar=None, degree=None, arg_real=None, arg_imag=None, *args, **kwargs):
    # Clenshaw summation of sinus with complex argument
    # [re, im] = clksin(ar, degree, arg_real, arg_imag);

    # Written by Kai Borre
    # December 20, 1995

    # See also WGS2UTM or CART2UTM

    # CVS record:
    # $Id: clksin.m,v 1.1.1.1.2.4 2006/08/22 13:45:59 dpl Exp $
    # ==========================================================================

    sin_arg_r = np.sin(arg_real)
    # clksin.m:14
    cos_arg_r = np.cos(arg_real)
    # clksin.m:15
    sinh_arg_i = np.sinh(arg_imag)
    # clksin.m:16
    cosh_arg_i = np.cosh(arg_imag)
    # clksin.m:17
    r = 2 * cos_arg_r * cosh_arg_i
    # clksin.m:19
    i = - 2 * sin_arg_r * sinh_arg_i
    # clksin.m:20
    hr1 = 0
    # clksin.m:22
    hr = 0
    # clksin.m:22
    hi1 = 0
    # clksin.m:22
    hi = 0
    # clksin.m:22
    # TODO fix index range of t
    for t in range(degree, 0, - 1):
        hr2 = hr1
        # clksin.m:25
        hr1 = hr
        # clksin.m:26
        hi2 = hi1
        # clksin.m:27
        hi1 = hi
        # clksin.m:28
        z = ar[t - 1] + r * hr1 - i * hi - hr2
        # clksin.m:29
        hi = i * hr1 + r * hi1 - hi2
        # clksin.m:30
        hr = z
    # clksin.m:31

    r = sin_arg_r * cosh_arg_i
    # clksin.m:34
    i = cos_arg_r * sinh_arg_i
    # clksin.m:35
    re = r * hr - i * hi
    # clksin.m:37
    im = r * hi + i * hr
    return re, im


# clksin.m:38

# cart2utm.m
def cart2utm(X=None, Y=None, Z=None, zone=None, *args, **kwargs):
    # CART2UTM  Transformation of (X,Y,Z) to (N,E,U) in UTM, zone 'zone'.

    # [E, N, U] = cart2utm(X, Y, Z, zone);

    #   Inputs:
    #       X,Y,Z       - Cartesian coordinates. Coordinates are referenced
    #                   with respect to the International Terrestrial Reference
    #                   Frame 1996 (ITRF96)
    #       zone        - UTM zone of the given position

    #   Outputs:
    #      E, N, U      - UTM coordinates (Easting, Northing, Uping)

    # Kai Borre -11-1994
    # Copyright (c) by Kai Borre

    # CVS record:
    # $Id: cart2utm.m,v 1.1.1.1.2.6 2007/01/30 09:45:12 dpl Exp $

    # This implementation is based upon
    # O. Andersson & K. Poder (1981) Koordinattransformationer
    #  ved Geod\ae{}tisk Institut. Landinspekt\oe{}ren
    #  Vol. 30: 552--571 and Vol. 31: 76

    # An excellent, general reference (KW) is
    # R. Koenig & K.H. Weise (1951) Mathematische Grundlagen der
    #  h\"oheren Geod\"asie und Kartographie.
    #  Erster Band, Springer Verlag

    # Explanation of variables used:
    # f	   flattening of ellipsoid
    # a	   semi major axis in m
    # m0	   1 - scale at central meridian; for UTM 0.0004
    # Q_n	   normalized meridian quadrant
    # E0	   Easting of central meridian
    # L0	   Longitude of central meridian
    # bg	   constants for ellipsoidal geogr. to spherical geogr.
    # gb	   constants for spherical geogr. to ellipsoidal geogr.
    # gtu	   constants for ellipsoidal N, E to spherical N, E
    # utg	   constants for spherical N, E to ellipoidal N, E
    # tolutm	tolerance for utm, 1.2E-10*meridian quadrant
    # tolgeo	tolerance for geographical, 0.00040 second of arc

    # B, L refer to latitude and longitude. Southern latitude is negative
    # International ellipsoid of 1924, valid for ED50

    a = 6378388.0
    # cart2utm.m:48
    f = 1.0 / 297.0
    # cart2utm.m:49
    ex2 = (2 - f) * f / (1 - f) ** 2
    # cart2utm.m:50
    c = a * np.sqrt(1 + ex2)
    # cart2utm.m:51
    vec = np.array([X, Y, Z - 4.5])
    # cart2utm.m:52
    alpha = 7.56e-07
    # cart2utm.m:53
    R = np.array([[1, - alpha, 0],
                  [alpha, 1, 0],
                  [0, 0, 1]])
    # cart2utm.m:54
    trans = np.array([89.5, 93.8, 127.6])
    # cart2utm.m:57
    scale = 0.9999988
    # cart2utm.m:58
    v = scale * R.dot(vec) + trans
    # cart2utm.m:59

    L = np.arctan2(v[1], v[0])
    # cart2utm.m:60
    N1 = 6395000.0
    # cart2utm.m:61

    B = np.arctan2(v[2] / ((1 - f) ** 2 * N1), np.linalg.norm(v[0:2]) / N1)
    # cart2utm.m:62

    U = 0.1
    # cart2utm.m:63
    oldU = 0
    # cart2utm.m:63
    iterations = 0
    # cart2utm.m:65
    while abs(U - oldU) > 0.0001:

        oldU = U
        # cart2utm.m:67
        N1 = c / np.sqrt(1 + ex2 * (np.cos(B)) ** 2)
        # cart2utm.m:68
        B = np.arctan2(v[2] / ((1 - f) ** 2 * N1 + U), np.linalg.norm(v[0:2]) / (N1 + U))
        # cart2utm.m:69
        U = np.linalg.norm(v[0:2]) / np.cos(B) - N1
        # cart2utm.m:70
        iterations += 1
        # cart2utm.m:72
        if iterations > 100:
            print 'Failed to approximate U with desired precision. U-oldU: %e.' % (U - oldU)
            break

    # Normalized meridian quadrant, KW p. 50 (96), p. 19 (38b), p. 5 (21)
    m0 = 0.0004
    # cart2utm.m:80
    n = f / (2 - f)
    # cart2utm.m:81
    m = n ** 2 * (1.0 / 4.0 + n ** 2 / 64)
    # cart2utm.m:82
    w = (a * (-n - m0 + m * (1 - m0))) / (1 + n)
    # cart2utm.m:83
    Q_n = a + w
    # cart2utm.m:84
    # Easting and longitude of central meridian
    E0 = 500000.0
    # cart2utm.m:87
    L0 = (zone - 30) * 6 - 3
    # cart2utm.m:88
    # Check tolerance for reverse transformation
    tolutm = np.pi / 2 * 1.2e-10 * Q_n
    # cart2utm.m:91
    tolgeo = 4e-05
    # cart2utm.m:92
    # Coefficients of trigonometric series

    # ellipsoidal to spherical geographical, KW p. 186--187, (51)-(52)
    # bg[1] = n*(-2 + n*(2/3    + n*(4/3	  + n*(-82/45))));
    # bg[2] = n^2*(5/3    + n*(-16/15 + n*(-13/9)));
    # bg[3] = n^3*(-26/15 + n*34/21);
    # bg[4] = n^4*1237/630;

    # spherical to ellipsoidal geographical, KW p. 190--191, (61)-(62)
    # gb[1] = n*(2	      + n*(-2/3    + n*(-2	 + n*116/45)));
    # gb[2] = n^2*(7/3    + n*(-8/5 + n*(-227/45)));
    # gb[3] = n^3*(56/15 + n*(-136/35));
    # gb[4] = n^4*4279/630;

    # spherical to ellipsoidal N, E, KW p. 196, (69)
    # gtu[1] = n*(1/2	  + n*(-2/3    + n*(5/16     + n*41/180)));
    # gtu[2] = n^2*(13/48	  + n*(-3/5 + n*557/1440));
    # gtu[3] = n^3*(61/240	 + n*(-103/140));
    # gtu[4] = n^4*49561/161280;

    # ellipsoidal to spherical N, E, KW p. 194, (65)
    # utg[1] = n*(-1/2	   + n*(2/3    + n*(-37/96	+ n*1/360)));
    # utg[2] = n^2*(-1/48	  + n*(-1/15 + n*437/1440));
    # utg[3] = n^3*(-17/480 + n*37/840);
    # utg[4] = n^4*(-4397/161280);

    # With f = 1/297 we get

    bg = np.array([- 0.00337077907, 4.73444769e-06, -8.2991457e-09, 1.5878533e-11])
    # cart2utm.m:122
    gb = np.array([0.00337077588, 6.6276908e-06, 1.78718601e-08, 5.49266312e-11])
    # cart2utm.m:127
    gtu = np.array([0.000841275991, 7.67306686e-07, 1.2129123e-09, 2.48508228e-12])
    # cart2utm.m:132
    utg = np.array([-0.000841276339, -5.95619298e-08, -1.69485209e-10, -2.20473896e-13])
    # cart2utm.m:137
    # Ellipsoidal latitude, longitude to spherical latitude, longitude
    neg_geo = False
    # cart2utm.m:143
    if B < 0:
        neg_geo = True
    # cart2utm.m:146

    Bg_r = np.abs(B)
    # cart2utm.m:149
    res_clensin = clsin(bg, 4, 2 * Bg_r)
    # cart2utm.m:150
    Bg_r = Bg_r + res_clensin
    # cart2utm.m:151
    L0 = L0 * np.pi / 180
    # cart2utm.m:152
    Lg_r = L - L0
    # cart2utm.m:153
    # Spherical latitude, longitude to complementary spherical latitude
    #  i.e. spherical N, E
    cos_BN = np.cos(Bg_r)
    # cart2utm.m:157
    Np = np.arctan2(np.sin(Bg_r), np.cos(Lg_r) * cos_BN)
    # cart2utm.m:158
    Ep = np.arctanh(np.sin(Lg_r) * cos_BN)
    # cart2utm.m:159
    # Spherical normalized N, E to ellipsoidal N, E
    Np *= 2
    # cart2utm.m:162
    Ep *= 2
    # cart2utm.m:163
    dN, dE = clksin(gtu, 4, Np, Ep)
    # cart2utm.m:164
    Np /= 2
    # cart2utm.m:165
    Ep /= 2
    # cart2utm.m:166
    Np += dN
    # cart2utm.m:167
    Ep += dE
    # cart2utm.m:168
    N = Q_n * Np
    # cart2utm.m:169
    E = Q_n * Ep + E0
    # cart2utm.m:170
    if neg_geo:
        N = -N + 20000000
    return E, N, U
# cart2utm.m:173

#################### end cart2utm.m ####################
