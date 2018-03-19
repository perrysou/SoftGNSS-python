import numpy as np


# check_t.m


def check_t(time=None, *args, **kwargs):
    # CHECK_T accounting for beginning or end of week crossover.

    # corrTime = check_t(time);

    #   Inputs:
    #       time        - time in seconds

    #   Outputs:
    #       corrTime    - corrected time (seconds)

    # Kai Borre 04-01-96
    # Copyright (c) by Kai Borre

    # CVS record:
    # $Id: check_t.m,v 1.1.1.1.2.4 2006/08/22 13:45:59 dpl Exp $
    # ==========================================================================

    half_week = 302400.0
    # check_t.m:19

    corrTime = time
    # check_t.m:21
    if time > half_week:
        corrTime = time - 2 * half_week
    # check_t.m:24
    elif time < - half_week:
        corrTime = time + 2 * half_week
    return corrTime


# check_t.m:26

####### end check_t.m  #################


# satpos.m


def satpos(transmitTime=None, prnList=None, eph=None, settings=None, *args, **kwargs):
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


    ## Initialize constants ===================================================
    numOfSatellites = prnList.size
    # satpos.m:29
    # GPS constatns

    gpsPi = 3.14159265359
    # satpos.m:33

    # system

    # --- Constants for satellite position calculation -------------------------
    Omegae_dot = 7.2921151467e-05
    # satpos.m:37

    GM = 3.986005e+14
    # satpos.m:38

    # the mass of the Earth, [m^3/s^2]
    F = - 4.442807633e-10
    # satpos.m:40

    ## Initialize results =====================================================
    satClkCorr = np.zeros(numOfSatellites)
    # satpos.m:43
    satPositions = np.zeros((3, numOfSatellites))
    # satpos.m:44
    ## Process each satellite =================================================

    for satNr in range(numOfSatellites):
        prn = prnList[satNr] - 1
        # satpos.m:50
        ## Find initial satellite clock correction --------------------------------
        # --- Find time difference ---------------------------------------------
        dt = check_t(transmitTime - eph[prn].t_oc)
        # satpos.m:55
        satClkCorr[satNr] = (eph[prn].a_f2 * dt + eph[prn].a_f1) * dt + eph[prn].a_f0 - eph[prn].T_GD
        # satpos.m:58
        time = transmitTime - satClkCorr[satNr]
        # satpos.m:62
        ## Find satellite's position ----------------------------------------------
        # Restore semi-major axis
        a = eph[prn].sqrtA * eph[prn].sqrtA
        # satpos.m:67
        tk = check_t(time - eph[prn].t_oe)
        # satpos.m:70
        n0 = np.sqrt(GM / a ** 3)
        # satpos.m:73
        n = n0 + eph[prn].deltan
        # satpos.m:75
        M = eph[prn].M_0 + n * tk
        # satpos.m:78
        M = np.remainder(M + 2 * gpsPi, 2 * gpsPi)
        # satpos.m:80
        E = M
        # satpos.m:83
        for ii in range(10):
            E_old = E
            # satpos.m:87
            E = M + eph[prn].e * np.sin(E)
            # satpos.m:88
            dE = np.remainder(E - E_old, 2 * gpsPi)
            # satpos.m:89
            if abs(dE) < 1e-12:
                # Necessary precision is reached, exit from the loop
                break
        # Reduce eccentric anomaly to between 0 and 360 deg
        E = np.remainder(E + 2 * gpsPi, 2 * gpsPi)
        # satpos.m:98
        dtr = F * eph[prn].e * eph[prn].sqrtA * np.sin(E)
        # satpos.m:101
        nu = np.arctan2(np.sqrt(1 - eph[prn].e ** 2) * np.sin(E), np.cos(E) - eph[prn].e)
        # satpos.m:104
        phi = nu + eph[prn].omega
        # satpos.m:107
        phi = np.remainder(phi, 2 * gpsPi)
        # satpos.m:109
        u = phi + eph[prn].C_uc * np.cos(2 * phi) + eph[prn].C_us * np.sin(2 * phi)
        # satpos.m:112
        r = a * (1 - eph[prn].e * np.cos(E)) + eph[prn].C_rc * np.cos(2 * phi) + eph[prn].C_rs * np.sin(2 * phi)
        # satpos.m:116
        i = eph[prn].i_0 + eph[prn].iDot * tk + eph[prn].C_ic * np.cos(2 * phi) + eph[prn].C_is * np.sin(2 * phi)
        # satpos.m:120
        Omega = eph[prn].omega_0 + (eph[prn].omegaDot - Omegae_dot) * tk - Omegae_dot * eph[prn].t_oe
        # satpos.m:125
        Omega = np.remainder(Omega + 2 * gpsPi, 2 * gpsPi)
        # satpos.m:128
        satPositions[0, satNr] = np.cos(u) * r * np.cos(Omega) - np.sin(u) * r * np.cos(i) * np.sin(Omega)
        # satpos.m:131
        satPositions[1, satNr] = np.cos(u) * r * np.sin(Omega) + np.sin(u) * r * np.cos(i) * np.cos(Omega)
        # satpos.m:132
        satPositions[2, satNr] = np.sin(u) * r * np.sin(i)
        # satpos.m:133
        ## Include relativistic correction in clock correction --------------------
        satClkCorr[satNr] = (eph[prn].a_f2 * dt + eph[prn].a_f1) * dt + eph[prn].a_f0 - eph[prn].T_GD + dtr
    return satPositions, satClkCorr
# satpos.m:137
