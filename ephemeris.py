def bin2dec(binaryStr):
    assert isinstance(binaryStr, str)
    return int(binaryStr, 2)


# twosComp2dec.m
def twosComp2dec(binaryStr):
    # TWOSCOMP2DEC(binaryNumber) Converts a two's-complement binary number
    # BINNUMBER (in Matlab it is a string type), represented as a row vector of
    # zeros and ones, to an integer.

    # intNumber = twosComp2dec(binaryNumber)

    # --- Check if the input is string -----------------------------------------
    if not isinstance(binaryStr, str):
        raise IOError('Input must be a string.')

    # --- Convert from binary form to a decimal number -------------------------
    intNumber = int(binaryStr, 2)

    # --- If the number was negative, then correct the result ------------------
    if binaryStr[0] == '1':
        intNumber -= 2 ** len(binaryStr)
    return intNumber


# checkPhase.m


def checkPhase(word, d30star):
    # Checks the parity of the supplied 30bit word.
    # The last parity bit of the previous word is used for the calculation.
    # A note on the procedure is supplied by the GPS standard positioning
    # service signal specification.

    # word = checkPhase(word, D30Star)

    #   Inputs:
    #       word        - an array with 30 bit long word from the navigation
    #                   message (a character array, must contain only '0' or
    #                   '1').
    #       D30Star     - the last bit of the previous word (char type).

    #   Outputs:
    #       word        - word with corrected polarity of the data bits
    #                   (character array).

    word_new = []
    if d30star == '1':
        # Data bits must be inverted
        for i in range(0, 24):
            if word[i] == '1':
                word[i] = '0'
            elif word[i] == '0':
                word[i] = '1'
    return word


# ephemeris.m
def ephemeris(bits, d30star):
    # Function decodes ephemerides and TOW from the given bit stream. The stream
    # (array) in the parameter BITS must contain 1500 bits. The first element in
    # the array must be the first bit of a subframe. The subframe ID of the
    # first subframe in the array is not important.

    # Function does not check parity!

    # [eph, TOW] = ephemeris(bits, D30Star)

    #   Inputs:
    #       bits        - bits of the navigation messages (5 subframes).
    #                   Type is character array and it must contain only
    #                   characters '0' or '1'.
    #       D30Star     - The last bit of the previous nav-word. Refer to the
    #                   GPS interface control document ICD (IS-GPS-200D) for
    #                   more details on the parity checking. Parameter type is
    #                   char. It must contain only characters '0' or '1'.
    #   Outputs:
    #       TOW         - Time Of Week (TOW) of the first sub-frame in the bit
    #                   stream (in seconds)
    #       eph         - SV ephemeris

    # Check if there is enough data ==========================================
    if len(bits) < 1500:
        raise TypeError('The parameter BITS must contain 1500 bits!')

    # Check if the parameters are strings ====================================
    if any([not isinstance(x, str) for x in bits]):
        raise TypeError('The parameter BITS must be a character array!')

    if not isinstance(d30star, str):
        raise TypeError('The parameter D30Star must be a char!')

    # Pi used in the GPS coordinate system
    gpsPi = 3.1415926535898

    # Decode all 5 sub-frames ================================================
    for i in range(5):
        # --- "Cut" one sub-frame's bits ---------------------------------------
        subframe = bits[300 * i:300 * (i + 1)]

        for j in range(10):
            subframe[30 * j: 30 * (j + 1)] = checkPhase(subframe[30 * j: 30 * (j + 1)], d30star)

            d30star = subframe[30 * (j + 1) - 1]

        # --- Decode the sub-frame id ------------------------------------------
        # For more details on sub-frame contents please refer to GPS IS.
        subframe = ''.join(subframe)
        subframeID = bin2dec(subframe[49:52])

        # The task is to select the necessary bits and convert them to decimal
        # numbers. For more details on sub-frame contents please refer to GPS
        # ICD (IS-GPS-200D).
        if 1 == subframeID:
            # It contains WN, SV clock corrections, health and accuracy
            weekNumber = bin2dec(subframe[60:70]) + 1024

            accuracy = bin2dec(subframe[72:76])

            health = bin2dec(subframe[76:82])

            T_GD = twosComp2dec(subframe[195:204]) * 2 ** (- 31)

            IODC = bin2dec(subframe[82:84] + subframe[196:204])

            t_oc = bin2dec(subframe[218:234]) * 2 ** 4

            a_f2 = twosComp2dec(subframe[240:248]) * 2 ** (- 55)

            a_f1 = twosComp2dec(subframe[248:264]) * 2 ** (- 43)

            a_f0 = twosComp2dec(subframe[270:292]) * 2 ** (- 31)

        elif 2 == subframeID:
            # It contains first part of ephemeris parameters
            IODE_sf2 = bin2dec(subframe[60:68])

            C_rs = twosComp2dec(subframe[68:84]) * 2 ** (- 5)

            deltan = twosComp2dec(subframe[90:106]) * 2 ** (- 43) * gpsPi

            M_0 = twosComp2dec(subframe[106:114] + subframe[120:144]) * 2 ** (- 31) * gpsPi

            C_uc = twosComp2dec(subframe[150:166]) * 2 ** (- 29)

            e = bin2dec(subframe[166:174] + subframe[180:204]) * 2 ** (- 33)

            C_us = twosComp2dec(subframe[210:226]) * 2 ** (- 29)

            sqrtA = bin2dec(subframe[226:234] + subframe[240:264]) * 2 ** (- 19)

            t_oe = bin2dec(subframe[270:286]) * 2 ** 4

        elif 3 == subframeID:
            # It contains second part of ephemeris parameters
            C_ic = twosComp2dec(subframe[60:76]) * 2 ** (- 29)

            omega_0 = twosComp2dec(subframe[76:84] + subframe[90:114]) * 2 ** (- 31) * gpsPi

            C_is = twosComp2dec(subframe[120:136]) * 2 ** (- 29)

            i_0 = twosComp2dec(subframe[136:144] + subframe[150:174]) * 2 ** (- 31) * gpsPi

            C_rc = twosComp2dec(subframe[180:196]) * 2 ** (- 5)

            omega = twosComp2dec(subframe[196:204] + subframe[210:234]) * 2 ** (- 31) * gpsPi

            omegaDot = twosComp2dec(subframe[240:264]) * 2 ** (- 43) * gpsPi

            IODE_sf3 = bin2dec(subframe[270:278])

            iDot = twosComp2dec(subframe[278:292]) * 2 ** (- 43) * gpsPi

        elif 4 == subframeID:
            # Almanac, ionospheric model, UTC parameters.
            # SV health (PRN: 25-32).
            # Not decoded at the moment.
            pass
        elif 5 == subframeID:
            # SV almanac and health (PRN: 1-24).
            # Almanac reference week number and time.
            # Not decoded at the moment.
            pass

    # Compute the time of week (TOW) of the first sub-frames in the array ====
    # Also correct the TOW. The transmitted TOW is actual TOW of the next
    # subframe and we need the TOW of the first subframe in this data block
    # (the variable subframe at this point contains bits of the last subframe).
    TOW = bin2dec(subframe[30:47]) * 6 - 30
    # Initialize fields for ephemeris
    eph = (weekNumber, accuracy, health, T_GD, IODC, t_oc, a_f2, a_f1, a_f0,
           IODE_sf2, C_rs, deltan, M_0, C_uc, e, C_us, sqrtA, t_oe,
           C_ic, omega_0, C_is, i_0, C_rc, omega, omegaDot, IODE_sf3, iDot)
    return eph, TOW
