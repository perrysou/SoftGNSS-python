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

    tropo = 0.0

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

        tropo = tropo + dr * ref * 1000

        if done:
            ddr = tropo

            break
        done = True

        refsea = (0.3719 / tksea - 1.292e-05) / tksea

        htop = 1.1385e-05 * (1255.0 / tksea + 0.05) / refsea

        ref = refsea * e0sea * ((htop - hsta) / htop) ** 4
    return ddr

######### end tropo.m  ###################
