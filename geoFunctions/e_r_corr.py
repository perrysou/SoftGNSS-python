import numpy as np


# e_r_corr.m


def e_r_corr(traveltime=None, X_sat=None, *args, **kwargs):
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

    # CVS record:
    # $Id: e_r_corr.m,v 1.1.1.1.2.6 2006/08/22 13:45:59 dpl Exp $
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
