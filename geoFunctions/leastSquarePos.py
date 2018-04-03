import numpy as np

from geoFunctions import e_r_corr, topocent, tropo


# leastSquarePos.m
def leastSquarePos(satpos=None, obs=None, settings=None, *args, **kwargs):
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

    X = satpos.copy()

    nmbOfSatellites = satpos.shape[1]

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

                Rot_X = e_r_corr.e_r_corr(traveltime, X[:, i])

                az[i], el[i], dist = topocent.topocent(pos[0:3], Rot_X - pos[0:3])

                if settings.useTropCorr:
                    # --- Calculate tropospheric correction --------------------
                    trop = tropo.tropo(np.sin(el[i] * dtr), 0.0, 1013.0, 293.0, 50.0, 0.0, 0.0, 0.0)

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
