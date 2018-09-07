import numpy as np


# cart2geo.m


def cart2geo(X, Y, Z, i, *args, **kwargs):
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

    # ==========================================================================

    a = np.array([6378388.0, 6378160.0, 6378135.0, 6378137.0, 6378137.0])

    f = np.array([1 / 297, 1 / 298.247, 1 / 298.26, 1 / 298.257222101, 1 / 298.257223563])

    lambda_ = np.arctan2(Y, X)

    ex2 = (2 - f[i]) * f[i] / ((1 - f[i]) ** 2)

    c = a[i] * np.sqrt(1 + ex2)

    phi = np.arctan(Z / (np.sqrt(X ** 2 + Y ** 2) * (1 - (2 - f[i])) * f[i]))

    h = 0.1

    oldh = 0

    iterations = 0

    while abs(h - oldh) > 1e-12:

        oldh = h

        N = c / np.sqrt(1 + ex2 * np.cos(phi) ** 2)

        phi = np.arctan(Z / (np.sqrt(X ** 2 + Y ** 2) * (1 - (2 - f[i]) * f[i] * N / (N + h))))

        h = np.sqrt(X ** 2 + Y ** 2) / np.cos(phi) - N

        iterations += 1

        if iterations > 100:
            print 'Failed to approximate h with desired precision. h-oldh: %e.' % (h - oldh)
            break

    phi *= (180 / np.pi)

    # b = zeros(1,3);
    # b(1,1) = fix(phi);
    # b(2,1) = fix(rem(phi,b(1,1))*60);
    # b(3,1) = (phi-b(1,1)-b(1,2)/60)*3600;

    lambda_ *= (180 / np.pi)

    # l = zeros(1,3);
    # l(1,1) = fix(lambda);
    # l(2,1) = fix(rem(lambda,l(1,1))*60);
    # l(3,1) = (lambda-l(1,1)-l(1,2)/60)*3600;

    # fprintf('\n     phi =#3.0f #3.0f #8.5f',b(1),b(2),b(3))
    # fprintf('\n  lambda =#3.0f #3.0f #8.5f',l(1),l(2),l(3))
    # fprintf('\n       h =#14.3f\n',h)
    return phi, lambda_, h


############## end cart2geo.m ###################


# clsin.m
def clsin(ar, degree, argument, *args, **kwargs):
    # Clenshaw summation of sinus of argument.

    # result = clsin(ar, degree, argument);

    # Written by Kai Borre
    # December 20, 1995

    # See also WGS2UTM or CART2UTM

    # ==========================================================================

    cos_arg = 2 * np.cos(argument)

    hr1 = 0

    hr = 0

    # TODO fix index range of t
    for t in range(degree, 0, -1):
        hr2 = hr1

        hr1 = hr

        hr = ar[t - 1] + cos_arg * hr1 - hr2

    result = hr * np.sin(argument)
    return result


####################### end clsin.m  #####################


# clksin.m
def clksin(ar, degree, arg_real, arg_imag, *args, **kwargs):
    # Clenshaw summation of sinus with complex argument
    # [re, im] = clksin(ar, degree, arg_real, arg_imag);

    # Written by Kai Borre
    # December 20, 1995

    # See also WGS2UTM or CART2UTM

    # ==========================================================================

    sin_arg_r = np.sin(arg_real)

    cos_arg_r = np.cos(arg_real)

    sinh_arg_i = np.sinh(arg_imag)

    cosh_arg_i = np.cosh(arg_imag)

    r = 2 * cos_arg_r * cosh_arg_i

    i = - 2 * sin_arg_r * sinh_arg_i

    hr1 = 0

    hr = 0

    hi1 = 0

    hi = 0

    # TODO fix index range of t
    for t in range(degree, 0, - 1):
        hr2 = hr1

        hr1 = hr

        hi2 = hi1

        hi1 = hi

        z = ar[t - 1] + r * hr1 - i * hi - hr2

        hi = i * hr1 + r * hi1 - hi2

        hr = z

    r = sin_arg_r * cosh_arg_i

    i = cos_arg_r * sinh_arg_i

    re = r * hr - i * hi

    im = r * hi + i * hr
    return re, im


# cart2utm.m
def cart2utm(X, Y, Z, zone, *args, **kwargs):
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

    f = 1.0 / 297.0

    ex2 = (2 - f) * f / (1 - f) ** 2

    c = a * np.sqrt(1 + ex2)

    vec = np.array([X, Y, Z - 4.5])

    alpha = 7.56e-07

    R = np.array([[1, - alpha, 0],
                  [alpha, 1, 0],
                  [0, 0, 1]])

    trans = np.array([89.5, 93.8, 127.6])

    scale = 0.9999988

    v = scale * R.dot(vec) + trans

    L = np.arctan2(v[1], v[0])

    N1 = 6395000.0

    B = np.arctan2(v[2] / ((1 - f) ** 2 * N1), np.linalg.norm(v[0:2]) / N1)

    U = 0.1

    oldU = 0

    iterations = 0

    while abs(U - oldU) > 0.0001:

        oldU = U

        N1 = c / np.sqrt(1 + ex2 * (np.cos(B)) ** 2)

        B = np.arctan2(v[2] / ((1 - f) ** 2 * N1 + U), np.linalg.norm(v[0:2]) / (N1 + U))

        U = np.linalg.norm(v[0:2]) / np.cos(B) - N1

        iterations += 1

        if iterations > 100:
            print 'Failed to approximate U with desired precision. U-oldU: %e.' % (U - oldU)
            break

    # Normalized meridian quadrant, KW p. 50 (96), p. 19 (38b), p. 5 (21)
    m0 = 0.0004

    n = f / (2 - f)

    m = n ** 2 * (1.0 / 4.0 + n ** 2 / 64)

    w = (a * (-n - m0 + m * (1 - m0))) / (1 + n)

    Q_n = a + w

    # Easting and longitude of central meridian
    E0 = 500000.0

    L0 = (zone - 30) * 6 - 3

    # Check tolerance for reverse transformation
    tolutm = np.pi / 2 * 1.2e-10 * Q_n

    tolgeo = 4e-05

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

    gb = np.array([0.00337077588, 6.6276908e-06, 1.78718601e-08, 5.49266312e-11])

    gtu = np.array([0.000841275991, 7.67306686e-07, 1.2129123e-09, 2.48508228e-12])

    utg = np.array([-0.000841276339, -5.95619298e-08, -1.69485209e-10, -2.20473896e-13])

    # Ellipsoidal latitude, longitude to spherical latitude, longitude
    neg_geo = False

    if B < 0:
        neg_geo = True

    Bg_r = np.abs(B)

    res_clensin = clsin(bg, 4, 2 * Bg_r)

    Bg_r = Bg_r + res_clensin

    L0 = L0 * np.pi / 180

    Lg_r = L - L0

    # Spherical latitude, longitude to complementary spherical latitude
    #  i.e. spherical N, E
    cos_BN = np.cos(Bg_r)

    Np = np.arctan2(np.sin(Bg_r), np.cos(Lg_r) * cos_BN)

    Ep = np.arctanh(np.sin(Lg_r) * cos_BN)

    # Spherical normalized N, E to ellipsoidal N, E
    Np *= 2

    Ep *= 2

    dN, dE = clksin(gtu, 4, Np, Ep)

    Np /= 2

    Ep /= 2

    Np += dN

    Ep += dE

    N = Q_n * Np

    E = Q_n * Ep + E0

    if neg_geo:
        N = -N + 20000000
    return E, N, U


#################### end cart2utm.m ####################


# deg2dms.m
def deg2dms(deg, *args, **kwargs):
    # DEG2DMS  Conversion of degrees to degrees, minutes, and seconds.
    # The output format (dms format) is: (degrees*100 + minutes + seconds/100)

    # Written by Kai Borre
    # February 7, 2001
    # Updated by Darius Plausinaitis

    # Save the sign for later processing
    neg_arg = False

    if deg < 0:
        # Only positive numbers should be used while spliting into deg/min/sec
        deg = -deg

        neg_arg = True

    # Split degrees minutes and seconds
    int_deg = np.floor(deg)

    decimal = deg - int_deg

    min_part = decimal * 60

    min_ = np.floor(min_part)

    sec_part = min_part - np.floor(min_part)

    sec = sec_part * 60

    # Check for overflow
    if sec == 60.0:
        min_ = min_ + 1

        sec = 0.0

    if min_ == 60.0:
        int_deg = int_deg + 1

        min_ = 0.0

    # Construct the output
    dmsOutput = int_deg * 100 + min_ + sec / 100

    # Correct the sign
    if neg_arg:
        dmsOutput = -dmsOutput
    return dmsOutput


################### end deg2dms.m ################


# dms2mat.m
def dms2mat(dmsInput, n, *args, **kwargs):
    # DMS2MAT  Splits a real a = dd*100 + mm + s/100 into[dd mm s.ssss]
    #         where n specifies the power of 10, to which the resulting seconds
    #         of the output should be rounded. E.g.: if a result is 23.823476
    #         seconds, and n = -3, then the output will be 23.823.

    # Written by Kai Borre
    # January 7, 2007
    # Updated by Darius Plausinaitis

    neg_arg = False

    if dmsInput < 0:
        # Only positive numbers should be used while spliting into deg/min/sec
        dmsInput = -dmsInput

        neg_arg = True

    # Split degrees minutes and seconds
    int_deg = np.floor(dmsInput / 100)

    mm = np.floor(dmsInput - 100 * int_deg)

    # we assume n<7; hence #2.10f is sufficient to hold ssdec
    ssdec = '%2.10f' % (dmsInput - 100 * int_deg - mm) * 100

    # Check for overflow
    if ssdec == 60.0:
        mm = mm + 1

        ssdec = 0.0

    if mm == 60.0:
        int_deg = int_deg + 1

        mm = 0.0

    # Corect the sign
    if neg_arg:
        int_deg = -int_deg

    # Compose the output
    matOutput = []
    matOutput[0] = int_deg

    matOutput[1] = mm

    matOutput[2] = float(ssdec[0:- n + 3])

    return matOutput


################### end dms2mat.m ################


# e_r_corr.m


def e_r_corr(traveltime, X_sat, *args, **kwargs):
    # E_R_CORR  Returns rotated satellite ECEF coordinates due to Earth
    # rotation during signal travel time

    # X_sat_rot = e_r_corr(traveltime, X_sat);

    #   Inputs:
    #       travelTime  - signal travel time
    #       X_sat       - satellite's ECEF coordinates

    #   Outputs:
    #       X_sat_rot   - rotated satellite's coordinates (ECEF)

    # Written by Kai Borre
    # Copyright (c) by Kai Borre

    # ==========================================================================

    Omegae_dot = 7.292115147e-05

    # --- Find rotation angle --------------------------------------------------
    omegatau = Omegae_dot * traveltime

    # --- Make a rotation matrix -----------------------------------------------
    R3 = np.array([[np.cos(omegatau), np.sin(omegatau), 0.0],
                   [-np.sin(omegatau), np.cos(omegatau), 0.0],
                   [0.0, 0.0, 1.0]])

    # --- Do the rotation ------------------------------------------------------
    X_sat_rot = R3.dot(X_sat)
    return X_sat_rot


######## end e_r_corr.m ####################

# findUtmZone.m


def findUtmZone(latitude, longitude, *args, **kwargs):
    # Function finds the UTM zone number for given longitude and latitude.
    # The longitude value must be between -180 (180 degree West) and 180 (180
    # degree East) degree. The latitude must be within -80 (80 degree South) and
    # 84 (84 degree North).

    # utmZone = findUtmZone(latitude, longitude);

    # Latitude and longitude must be in decimal degrees (e.g. 15.5 degrees not
    # 15 deg 30 min).

    # Check value bounds =====================================================

    if longitude > 180 or longitude < - 180:
        raise IOError('Longitude value exceeds limits (-180:180).')

    if latitude > 84 or latitude < - 80:
        raise IOError('Latitude value exceeds limits (-80:84).')

    # Find zone ==============================================================

    # Start at 180 deg west = -180 deg

    utmZone = np.fix((180 + longitude) / 6) + 1

    # Correct zone numbers for particular areas ==============================

    if latitude > 72:
        # Corrections for zones 31 33 35 37
        if 0 <= longitude < 9:
            utmZone = 31

        elif 9 <= longitude < 21:
            utmZone = 33

        elif 21 <= longitude < 33:
            utmZone = 35

        elif 33 <= longitude < 42:
            utmZone = 37

    elif 56 <= latitude < 64:
        # Correction for zone 32
        if 3 <= longitude < 12:
            utmZone = 32
    return utmZone


# geo2cart.m
def geo2cart(phi, lambda_, h, i=4, *args, **kwargs):
    # GEO2CART Conversion of geographical coordinates (phi, lambda, h) to
    # Cartesian coordinates (X, Y, Z).

    # [X, Y, Z] = geo2cart(phi, lambda, h, i);

    # Format for phi and lambda: [degrees minutes seconds].
    # h, X, Y, and Z are in meters.

    # Choices i of Reference Ellipsoid
    #   1. International Ellipsoid 1924
    #   2. International Ellipsoid 1967
    #   3. World Geodetic System 1972
    #   4. Geodetic Reference System 1980
    #   5. World Geodetic System 1984

    #   Inputs:
    #       phi       - geocentric latitude (format [degrees minutes seconds])
    #       lambda    - geocentric longitude (format [degrees minutes seconds])
    #       h         - height
    #       i         - reference ellipsoid type

    #   Outputs:
    #       X, Y, Z   - Cartesian coordinates (meters)

    # Kai Borre 10-13-98
    # Copyright (c) by Kai Borre

    # ==========================================================================

    b = phi[0] + phi[1] / 60. + phi[2] / 3600.

    b = b * np.pi / 180.

    l = lambda_[1] + lambda_[2] / 60. + lambda_[3] / 3600.

    l = l * np.pi / 180.

    a = [6378388, 6378160, 6378135, 6378137, 6378137]

    f = [1 / 297, 1 / 298.247, 1 / 298.26, 1 / 298.257222101, 1 / 298.257223563]

    ex2 = (2. - f[i]) * f[i] / (1. - f[i]) ** 2

    c = a[i] * np.sqrt(1. + ex2)

    N = c / np.sqrt(1. + ex2 * np.cos(b) ** 2)

    X = (N + h) * np.cos(b) * np.cos(l)

    Y = (N + h) * np.cos(b) * np.sin(l)

    Z = ((1. - f[i]) ** 2 * N + h) * np.sin(b)

    return X, Y, Z


# leastSquarePos.m
def leastSquarePos(satpos_, obs, settings, *args, **kwargs):
    # Function calculates the Least Square Solution.

    # [pos, el, az, dop] = leastSquarePos(satpos, obs, settings);

    #   Inputs:
    #       satpos      - Satellites positions (in ECEF system: [X; Y; Z;] -
    #                   one column per satellite)
    #       obs         - Observations - the pseudorange measurements to each
    #                   satellite:
    #                   (e.g. [20000000 21000000 .... .... .... .... ....])
    #       settings    - receiver settings

    #   Outputs:
    #       pos         - receiver position and receiver clock error
    #                   (in ECEF system: [X, Y, Z, dt])
    #       el          - Satellites elevation angles (degrees)
    #       az          - Satellites azimuth angles (degrees)
    #       dop         - Dilutions Of Precision ([GDOP PDOP HDOP VDOP TDOP])

    # === Initialization =======================================================
    nmbOfIterations = 7

    dtr = np.pi / 180

    pos = np.zeros(4)

    X = satpos_.copy()

    nmbOfSatellites = satpos_.shape[1]

    A = np.zeros((nmbOfSatellites, 4))

    omc = np.zeros(nmbOfSatellites)

    az = np.zeros(nmbOfSatellites)

    el = np.zeros(nmbOfSatellites)

    dop = np.zeros(5)
    # === Iteratively find receiver position ===================================
    for iter_ in range(nmbOfIterations):
        for i in range(nmbOfSatellites):
            if iter_ == 0:
                # --- Initialize variables at the first iteration --------------
                Rot_X = X[:, i].copy()

                trop = 2

            else:
                # --- Update equations -----------------------------------------
                rho2 = (X[0, i] - pos[0]) ** 2 + (X[1, i] - pos[1]) ** 2 + (X[2, i] - pos[2]) ** 2

                traveltime = np.sqrt(rho2) / settings.c

                Rot_X = e_r_corr(traveltime, X[:, i])

                az[i], el[i], dist = topocent(pos[0:3], Rot_X - pos[0:3])

                if settings.useTropCorr:
                    # --- Calculate tropospheric correction --------------------
                    trop = tropo(np.sin(el[i] * dtr), 0.0, 1013.0, 293.0, 50.0, 0.0, 0.0, 0.0)

                else:
                    # Do not calculate or apply the tropospheric corrections
                    trop = 0

            # --- Apply the corrections ----------------------------------------
            omc[i] = obs[i] - np.linalg.norm(Rot_X - pos[0:3]) - pos[3] - trop

            A[i, :] = np.array([-(Rot_X[0] - pos[0]) / obs[i],
                                -(Rot_X[1] - pos[1]) / obs[i],
                                -(Rot_X[2] - pos[2]) / obs[i],
                                1])

        # These lines allow the code to exit gracefully in case of any errors
        if np.linalg.matrix_rank(A) != 4:
            pos = np.zeros((4, 1))

            return pos, el, az, dop
        # --- Find position update ---------------------------------------------
        x = np.linalg.lstsq(A, omc, rcond=None)[0]

        pos = pos + x.flatten()

    pos = pos.T

    # === Calculate Dilution Of Precision ======================================
    # --- Initialize output ------------------------------------------------
    # dop = np.zeros((1, 5))

    Q = np.linalg.inv(A.T.dot(A))

    dop[0] = np.sqrt(np.trace(Q))

    dop[1] = np.sqrt(Q[0, 0] + Q[1, 1] + Q[2, 2])

    dop[2] = np.sqrt(Q[0, 0] + Q[1, 1])

    dop[3] = np.sqrt(Q[2, 2])

    dop[4] = np.sqrt(Q[3, 3])

    return pos, el, az, dop


# check_t.m


def check_t(time, *args, **kwargs):
    # CHECK_T accounting for beginning or end of week crossover.

    # corrTime = check_t(time);

    #   Inputs:
    #       time        - time in seconds

    #   Outputs:
    #       corrTime    - corrected time (seconds)

    # Kai Borre 04-01-96
    # Copyright (c) by Kai Borre

    # ==========================================================================

    half_week = 302400.0

    corrTime = time

    if time > half_week:
        corrTime = time - 2 * half_week

    elif time < - half_week:
        corrTime = time + 2 * half_week
    return corrTime


####### end check_t.m  #################


# satpos.m


def satpos(transmitTime, prnList, eph, settings, *args, **kwargs):
    # SATPOS Computation of satellite coordinates X,Y,Z at TRANSMITTIME for
    # given ephemeris EPH. Coordinates are computed for each satellite in the
    # list PRNLIST.
    # [satPositions, satClkCorr] = satpos(transmitTime, prnList, eph, settings);

    #   Inputs:
    #       transmitTime  - transmission time
    #       prnList       - list of PRN-s to be processed
    #       eph           - ephemerides of satellites
    #       settings      - receiver settings

    #   Outputs:
    #       satPositions  - position of satellites (in ECEF system [X; Y; Z;])
    #       satClkCorr    - correction of satellite clocks

    # Initialize constants ===================================================
    numOfSatellites = prnList.size

    # GPS constatns

    gpsPi = 3.14159265359

    # system

    # --- Constants for satellite position calculation -------------------------
    Omegae_dot = 7.2921151467e-05

    GM = 3.986005e+14

    # the mass of the Earth, [m^3/s^2]
    F = - 4.442807633e-10

    # Initialize results =====================================================
    satClkCorr = np.zeros(numOfSatellites)

    satPositions = np.zeros((3, numOfSatellites))

    # Process each satellite =================================================

    for satNr in range(numOfSatellites):
        prn = prnList[satNr] - 1

        # Find initial satellite clock correction --------------------------------
        # --- Find time difference ---------------------------------------------
        dt = check_t(transmitTime - eph[prn].t_oc)

        satClkCorr[satNr] = (eph[prn].a_f2 * dt + eph[prn].a_f1) * dt + eph[prn].a_f0 - eph[prn].T_GD

        time = transmitTime - satClkCorr[satNr]

        # Find satellite's position ----------------------------------------------
        # Restore semi-major axis
        a = eph[prn].sqrtA * eph[prn].sqrtA

        tk = check_t(time - eph[prn].t_oe)

        n0 = np.sqrt(GM / a ** 3)

        n = n0 + eph[prn].deltan

        M = eph[prn].M_0 + n * tk

        M = np.remainder(M + 2 * gpsPi, 2 * gpsPi)

        E = M

        for ii in range(10):
            E_old = E

            E = M + eph[prn].e * np.sin(E)

            dE = np.remainder(E - E_old, 2 * gpsPi)

            if abs(dE) < 1e-12:
                # Necessary precision is reached, exit from the loop
                break
        # Reduce eccentric anomaly to between 0 and 360 deg
        E = np.remainder(E + 2 * gpsPi, 2 * gpsPi)

        dtr = F * eph[prn].e * eph[prn].sqrtA * np.sin(E)

        nu = np.arctan2(np.sqrt(1 - eph[prn].e ** 2) * np.sin(E), np.cos(E) - eph[prn].e)

        phi = nu + eph[prn].omega

        phi = np.remainder(phi, 2 * gpsPi)

        u = phi + eph[prn].C_uc * np.cos(2 * phi) + eph[prn].C_us * np.sin(2 * phi)

        r = a * (1 - eph[prn].e * np.cos(E)) + eph[prn].C_rc * np.cos(2 * phi) + eph[prn].C_rs * np.sin(2 * phi)

        i = eph[prn].i_0 + eph[prn].iDot * tk + eph[prn].C_ic * np.cos(2 * phi) + eph[prn].C_is * np.sin(2 * phi)

        Omega = eph[prn].omega_0 + (eph[prn].omegaDot - Omegae_dot) * tk - Omegae_dot * eph[prn].t_oe

        Omega = np.remainder(Omega + 2 * gpsPi, 2 * gpsPi)

        satPositions[0, satNr] = np.cos(u) * r * np.cos(Omega) - np.sin(u) * r * np.cos(i) * np.sin(Omega)

        satPositions[1, satNr] = np.cos(u) * r * np.sin(Omega) + np.sin(u) * r * np.cos(i) * np.cos(Omega)

        satPositions[2, satNr] = np.sin(u) * r * np.sin(i)

        # Include relativistic correction in clock correction --------------------
        satClkCorr[satNr] = (eph[prn].a_f2 * dt + eph[prn].a_f1) * dt + eph[prn].a_f0 - eph[prn].T_GD + dtr
    return satPositions, satClkCorr


# topocent.m


# togeod.m
def togeod(a, finv, X, Y, Z, *args, **kwargs):
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


def topocent(X, dx, *args, **kwargs):
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


# tropo.m


def tropo(sinel, hsta, p, tkel, hum, hp, htkel, hhum):
    # TROPO  Calculation of tropospheric correction.
    #       The range correction ddr in m is to be subtracted from
    #       pseudo-ranges and carrier phases

    # ddr = tropo(sinel, hsta, p, tkel, hum, hp, htkel, hhum);

    #   Inputs:
    #       sinel   - sin of elevation angle of satellite
    #       hsta    - height of station in km
    #       p       - atmospheric pressure in mb at height hp
    #       tkel    - surface temperature in degrees Kelvin at height htkel
    #       hum     - humidity in # at height hhum
    #       hp      - height of pressure measurement in km
    #       htkel   - height of temperature measurement in km
    #       hhum    - height of humidity measurement in km

    #   Outputs:
    #       ddr     - range correction (meters)

    # Reference
    # Goad, C.C. & Goodman, L. (1974) A Modified Tropospheric
    # Refraction Correction Model. Paper presented at the
    # American Geophysical Union Annual Fall Meeting, San
    # Francisco, December 12-17

    # A Matlab reimplementation of a C code from driver.
    # Kai Borre 06-28-95

    # ==========================================================================

    a_e = 6378.137

    b0 = 7.839257e-05

    tlapse = -6.5

    tkhum = tkel + tlapse * (hhum - htkel)

    atkel = 7.5 * (tkhum - 273.15) / (237.3 + tkhum - 273.15)

    e0 = 0.0611 * hum * 10 ** atkel

    tksea = tkel - tlapse * htkel

    em = -978.77 / (2870400.0 * tlapse * 1e-05)

    tkelh = tksea + tlapse * hhum

    e0sea = e0 * (tksea / tkelh) ** (4 * em)

    tkelp = tksea + tlapse * hp

    psea = p * (tksea / tkelp) ** em

    if sinel < 0:
        sinel = 0

    tropo_ = 0.0

    done = False

    refsea = 7.7624e-05 / tksea

    htop = 1.1385e-05 / refsea

    refsea = refsea * psea

    ref = refsea * ((htop - hsta) / htop) ** 4

    while 1:
        rtop = (a_e + htop) ** 2 - (a_e + hsta) ** 2 * (1 - sinel ** 2)

        if rtop < 0:
            rtop = 0

        rtop = np.sqrt(rtop) - (a_e + hsta) * sinel

        a = -sinel / (htop - hsta)

        b = -b0 * (1 - sinel ** 2) / (htop - hsta)

        rn = np.zeros(8)

        for i in range(8):
            rn[i] = rtop ** (i + 2)

        alpha = np.array([2 * a, 2 * a ** 2 + 4 * b / 3, a * (a ** 2 + 3 * b),
                          a ** 4 / 5 + 2.4 * a ** 2 * b + 1.2 * b ** 2,
                          2 * a * b * (a ** 2 + 3 * b) / 3,
                          b ** 2 * (6 * a ** 2 + 4 * b) * 0.1428571, 0, 0])

        if b ** 2 > 1e-35:
            alpha[6] = a * b ** 3 / 2

            alpha[7] = b ** 4 / 9

        dr = rtop

        dr = dr + alpha.dot(rn)

        tropo_ += dr * ref * 1000

        if done:
            ddr = tropo_

            break
        done = True

        refsea = (0.3719 / tksea - 1.292e-05) / tksea

        htop = 1.1385e-05 * (1255.0 / tksea + 0.05) / refsea

        ref = refsea * e0sea * ((htop - hsta) / htop) ** 4
    return ddr


######### end tropo.m  ###################


if __name__ == '__main__':
    # print "This program is being run by itself"
    pass
else:
    print 'Importing functions from ./geoFunctions/'
