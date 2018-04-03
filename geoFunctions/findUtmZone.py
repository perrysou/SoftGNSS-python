# findUtmZone.m
import numpy as np


def findUtmZone(latitude=None, longitude=None, *args, **kwargs):
    # Function finds the UTM zone number for given longitude and latitude.
    # The longitude value must be between -180 (180 degree West) and 180 (180
    # degree East) degree. The latitude must be within -80 (80 degree South) and
    # 84 (84 degree North).

    # utmZone = findUtmZone(latitude, longitude);

    # Latitude and longitude must be in decimal degrees (e.g. 15.5 degrees not
    # 15 deg 30 min).

    ## Check value bounds =====================================================

    if longitude > 180 or longitude < - 180:
        raise IOError('Longitude value exceeds limits (-180:180).')

    if latitude > 84 or latitude < - 80:
        raise IOError('Latitude value exceeds limits (-80:84).')

    ## Find zone ==============================================================

    # Start at 180 deg west = -180 deg

    utmZone = np.fix((180 + longitude) / 6) + 1

    ## Correct zone numbers for particular areas ==============================

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
