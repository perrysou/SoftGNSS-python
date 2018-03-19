import numpy as np


# tropo.m


def tropo(sinel=None, hsta=None, p=None, tkel=None, hum=None, hp=None, htkel=None, hhum=None, *args, **kwargs):
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

    # CVS record:
    # $Id: tropo.m,v 1.1.1.1.2.4 2006/08/22 13:46:00 dpl Exp $
    # ==========================================================================

    a_e = 6378.137
    # tropo.m:34

    b0 = 7.839257e-05
    # tropo.m:35
    tlapse = -6.5
    # tropo.m:36
    tkhum = tkel + tlapse * (hhum - htkel)
    # tropo.m:37
    atkel = 7.5 * (tkhum - 273.15) / (237.3 + tkhum - 273.15)
    # tropo.m:38
    e0 = 0.0611 * hum * 10 ** atkel
    # tropo.m:39
    tksea = tkel - tlapse * htkel
    # tropo.m:40
    em = -978.77 / (2870400.0 * tlapse * 1e-05)
    # tropo.m:41
    tkelh = tksea + tlapse * hhum
    # tropo.m:42
    e0sea = e0 * (tksea / tkelh) ** (4 * em)
    # tropo.m:43
    tkelp = tksea + tlapse * hp
    # tropo.m:44
    psea = p * (tksea / tkelp) ** em
    # tropo.m:45
    if sinel < 0:
        sinel = 0
    # tropo.m:48

    tropo = 0.0
    # tropo.m:51
    done = False
    # tropo.m:52
    refsea = 7.7624e-05 / tksea
    # tropo.m:53
    htop = 1.1385e-05 / refsea
    # tropo.m:54
    refsea = refsea * psea
    # tropo.m:55
    ref = refsea * ((htop - hsta) / htop) ** 4
    # tropo.m:56
    while 1:
        rtop = (a_e + htop) ** 2 - (a_e + hsta) ** 2 * (1 - sinel ** 2)
        # tropo.m:59
        if rtop < 0:
            rtop = 0
        # tropo.m:63
        rtop = np.sqrt(rtop) - (a_e + hsta) * sinel
        # tropo.m:66
        a = -sinel / (htop - hsta)
        # tropo.m:67
        b = -b0 * (1 - sinel ** 2) / (htop - hsta)
        # tropo.m:68
        rn = np.zeros(8)
        # tropo.m:69
        for i in range(8):
            rn[i] = rtop ** (i + 2)
        # tropo.m:72
        alpha = np.array([2 * a, 2 * a ** 2 + 4 * b / 3, a * (a ** 2 + 3 * b),
                           a ** 4 / 5 + 2.4 * a ** 2 * b + 1.2 * b ** 2,
                           2 * a * b * (a ** 2 + 3 * b) / 3,
                           b ** 2 * (6 * a ** 2 + 4 * b) * 0.1428571, 0, 0])
        # tropo.m:75
        if b ** 2 > 1e-35:
            alpha[6] = a * b ** 3 / 2
            # tropo.m:80
            alpha[7] = b ** 4 / 9
        # tropo.m:81
        dr = rtop
        # tropo.m:84
        dr = dr + alpha.dot(rn)
        # tropo.m:85
        tropo = tropo + dr * ref * 1000
        # tropo.m:86
        if done:
            ddr = tropo
            # tropo.m:89
            break
        done = True
        # tropo.m:93
        refsea = (0.3719 / tksea - 1.292e-05) / tksea
        # tropo.m:94
        htop = 1.1385e-05 * (1255.0 / tksea + 0.05) / refsea
        # tropo.m:95
        ref = refsea * e0sea * ((htop - hsta) / htop) ** 4
    return ddr
# tropo.m:96


######### end tropo.m  ###################
