import numpy as np

import ephemeris
import findPreambles
import initSettings
from geoFunctions import satpos, leastSquarePos, cart2geo, findUtmZone, cart2utm


# ./calculatePseudoranges.m
def calculatePseudoranges(trackResults=None, msOfTheSignal=None, channelList=None, settings=None, *args, **kwargs):
    # calculatePseudoranges finds relative pseudoranges for all satellites
    # listed in CHANNELLIST at the specified millisecond of the processed
    # signal. The pseudoranges contain unknown receiver clock offset. It can be
    # found by the least squares position search procedure.

    # [pseudoranges] = calculatePseudoranges(trackResults, msOfTheSignal, ...
    #                                       channelList, settings)

    #   Inputs:
    #       trackResults    - output from the tracking function
    #       msOfTheSignal   - pseudorange measurement point (millisecond) in
    #                       the trackResults structure
    #       channelList     - list of channels to be processed
    #       settings        - receiver settings

    #   Outputs:
    #       pseudoranges    - relative pseudoranges to the satellites.

    # --- Set initial travel time to infinity ----------------------------------
    # Later in the code a shortest pseudorange will be selected. Therefore
    # pseudoranges from non-tracking channels must be the longest - e.g.
    # infinite.
    travelTime = np.Inf * np.ones(settings.numberOfChannels)

    # Find number of samples per spreading code
    samplesPerCode = np.round(settings.samplingFreq / (settings.codeFreqBasis / settings.codeLength))

    # --- For all channels in the list ...
    for channelNr in channelList:
        # --- Compute the travel times -----------------------------------------
        travelTime[channelNr] = trackResults[channelNr].absoluteSample[
                                    np.int(msOfTheSignal[channelNr])] / samplesPerCode

    # --- Truncate the travelTime and compute pseudoranges ---------------------
    minimum = np.floor(travelTime.min())

    travelTime = travelTime - minimum + settings.startOffset

    # --- Convert travel time to a distance ------------------------------------
    # The speed of light must be converted from meters per second to meters
    # per millisecond.
    pseudoranges = travelTime * settings.c / 1000
    return pseudoranges


# ./postNavigation.m
def postNavigation(trackResults=None, settings=None, *args, **kwargs):
    # Function calculates navigation solutions for the receiver (pseudoranges,
    # positions). At the end it converts coordinates from the WGS84 system to
    # the UTM, geocentric or any additional coordinate system.

    # [navSolutions, eph] = postNavigation(trackResults, settings)

    #   Inputs:
    #       trackResults    - results from the tracking function (structure
    #                       array).
    #       settings        - receiver settings.
    #   Outputs:
    #       navSolutions    - contains measured pseudoranges, receiver
    #                       clock error, receiver coordinates in several
    #                       coordinate systems (at least ECEF and UTM).
    #       eph             - received ephemerides of all SV (structure array).

    ## Check is there enough data to obtain any navigation solution ===========
    # It is necessary to have at least three subframes (number 1, 2 and 3) to
    # find satellite coordinates. Then receiver position can be found too.
    # The function requires all 5 subframes, because the tracking starts at
    # arbitrary point. Therefore the first received subframes can be any three
    # from the 5.
    # One subframe length is 6 seconds, therefore we need at least 30 sec long
    # record (5 * 6 = 30 sec = 30000ms). We add extra seconds for the cases,
    # when tracking has started in a middle of a subframe.

    if settings.msToProcess < 36000 or sum(trackResults.status != '-') < 4:
        # Show the error message and exit
        print 'Record is to short or too few satellites tracked. Exiting!'
        navSolutions = []

        eph = []

        return navSolutions, eph

    ## Find preamble start positions ==========================================

    subFrameStart, activeChnList = findPreambles.findPreambles(trackResults, settings)

    ## Decode ephemerides =====================================================
    eph = np.recarray((32,), formats=['O'] * 27,
                      names='weekNumber,accuracy,health,T_GD,IODC,t_oc,a_f2,a_f1,a_f0,IODE_sf2,C_rs,deltan,M_0,C_uc,e,C_us,sqrtA,t_oe,C_ic,omega_0,C_is,i_0,C_rc,omega,omegaDot,IODE_sf3,iDot')
    for channelNr in activeChnList:
        # === Convert tracking output to navigation bits =======================
        # --- Copy 5 sub-frames long record from tracking output ---------------
        navBitsSamples = trackResults[channelNr].I_P[subFrameStart[channelNr] - 20:
                                                     subFrameStart[channelNr] + 1500 * 20].copy()

        navBitsSamples = navBitsSamples.reshape(20, -1, order='F')

        navBits = navBitsSamples.sum(0)

        # The expression (navBits > 0) returns an array with elements set to 1
        # if the condition is met and set to 0 if it is not met.
        navBits = (navBits > 0) * 1

        # The function ephemeris expects input in binary form. In Matlab it is
        # a string array containing only "0" and "1" characters.
        navBitsBin = map(str, navBits)

        eph[trackResults[channelNr].PRN - 1], TOW = ephemeris.ephemeris(navBitsBin[1:], navBitsBin[0])

        if eph[trackResults[channelNr].PRN - 1].IODC is None or \
                eph[trackResults[channelNr].PRN - 1].IODE_sf2 is None or \
                eph[trackResults[channelNr].PRN - 1].IODE_sf3 is None:
            # --- Exclude channel from the list (from further processing) ------
            activeChnList = np.setdiff1d(activeChnList, channelNr)

    ## Check if the number of satellites is still above 3 =====================
    if activeChnList.size == 0 or activeChnList.size < 4:
        # Show error message and exit
        print 'Too few satellites with ephemeris data for postion calculations. Exiting!'
        navSolutions = None

        eph = None

        return navSolutions, eph

    ## Initialization =========================================================

    # Set the satellite elevations array to INF to include all satellites for
    # the first calculation of receiver position. There is no reference point
    # to find the elevation angle as there is no receiver position estimate at
    # this point.
    satElev = np.Inf * np.ones(settings.numberOfChannels)

    # Save the active channel list. The list contains satellites that are
    # tracked and have the required ephemeris data. In the next step the list
    # will depend on each satellite's elevation angle, which will change over
    # time.
    readyChnList = activeChnList

    transmitTime = TOW

    ###########################################################################
    ##   Do the satellite and receiver position calculations                  #
    ###########################################################################
    ## Initialization of current measurement ==================================
    channel = np.rec.array([(np.zeros((settings.numberOfChannels, 64)),
                             np.nan * np.ones((settings.numberOfChannels, 64)),
                             np.nan * np.ones((settings.numberOfChannels, 64)),
                             np.nan * np.ones((settings.numberOfChannels, 64)),
                             np.nan * np.ones((settings.numberOfChannels, 64))
                             )], formats=['O'] * 5, names='PRN,el,az,rawP,correctedP')
    navSolutions = np.rec.array([(channel,
                                  np.zeros((5, 64)),
                                  np.nan * np.ones(64),
                                  np.nan * np.ones(64),
                                  np.nan * np.ones(64),
                                  np.nan * np.ones(64),
                                  np.nan * np.ones(64),
                                  np.nan * np.ones(64),
                                  np.nan * np.ones(64),
                                  0,
                                  np.nan * np.ones(64),
                                  np.nan * np.ones(64),
                                  np.nan * np.ones(64)
                                  )], formats=['O'] * 13,
                                names='channel,DOP,X,Y,Z,dt,latitude,longitude,height,utmZone,E,N,U')
    for currMeasNr in range(np.int(np.fix(settings.msToProcess - subFrameStart.max()) / settings.navSolPeriod)):
        # Exclude satellites, that are belove elevation mask
        activeChnList = np.intersect1d((satElev >= settings.elevationMask).nonzero()[0], readyChnList)

        channel[0].PRN[activeChnList, currMeasNr] = trackResults[activeChnList].PRN

        # do to elevation mask will not "jump" to position (0,0) in the sky
        # plot.
        # channel[0].el[:, currMeasNr] = np.nan * np.ones(settings.numberOfChannels)

        # channel[0].az[:, currMeasNr] = np.nan * np.ones(settings.numberOfChannels)

        ## Find pseudoranges ======================================================
        channel[0].rawP[:, currMeasNr] = calculatePseudoranges(trackResults,
                                                               subFrameStart + settings.navSolPeriod * currMeasNr,
                                                               activeChnList, settings)

        ## Find satellites positions and clocks corrections =======================
        satPositions, satClkCorr = satpos.satpos(transmitTime, trackResults[activeChnList].PRN, eph, settings)

        ## Find receiver position =================================================
        # 3D receiver position can be found only if signals from more than 3
        # satellites are available
        if activeChnList.size > 3:
            # === Calculate receiver position ==================================
            (xyzdt,
             channel[0].el[activeChnList, currMeasNr],
             channel[0].az[activeChnList, currMeasNr],
             navSolutions[0].DOP[:, currMeasNr]) = leastSquarePos.leastSquarePos(satPositions,
                                                                                 channel[0].rawP[
                                                                                     activeChnList, currMeasNr] +
                                                                                 satClkCorr * settings.c,
                                                                                 settings)

            navSolutions[0].X[currMeasNr] = xyzdt[0]

            navSolutions[0].Y[currMeasNr] = xyzdt[1]

            navSolutions[0].Z[currMeasNr] = xyzdt[2]

            navSolutions[0].dt[currMeasNr] = xyzdt[3]

            satElev = channel[0].el[:, currMeasNr]

            channel[0].correctedP[activeChnList, currMeasNr] = channel[0].rawP[activeChnList, currMeasNr] + \
                                                               satClkCorr * settings.c + \
                                                               navSolutions[0].dt[currMeasNr]

            ## Coordinate conversion ==================================================
            # === Convert to geodetic coordinates ==============================
            (navSolutions[0].latitude[currMeasNr],
             navSolutions[0].longitude[currMeasNr],
             navSolutions[0].height[currMeasNr]) = cart2geo.cart2geo(navSolutions[0].X[currMeasNr],
                                                                     navSolutions[0].Y[currMeasNr],
                                                                     navSolutions[0].Z[currMeasNr],
                                                                     4)

            navSolutions[0].utmZone = findUtmZone.findUtmZone(navSolutions[0].latitude[currMeasNr],
                                                              navSolutions[0].longitude[currMeasNr])

            (navSolutions[0].E[currMeasNr],
             navSolutions[0].N[currMeasNr],
             navSolutions[0].U[currMeasNr]) = cart2utm.cart2utm(xyzdt[0], xyzdt[1], xyzdt[2],
                                                                navSolutions[0].utmZone)

        else:
            # --- There are not enough satellites to find 3D position ----------
            print '   Measurement No. %d' % currMeasNr + ': Not enough information for position solution.'
            # excluded automatically in all plots. For DOP it is easier to use
            # zeros. NaN values might need to be excluded from results in some
            # of further processing to obtain correct results.
            navSolutions[0].X[currMeasNr] = np.nan

            navSolutions[0].Y[currMeasNr] = np.nan

            navSolutions[0].Z[currMeasNr] = np.nan

            navSolutions[0].dt[currMeasNr] = np.nan

            navSolutions[0].DOP[:, currMeasNr] = np.zeros(5)

            navSolutions[0].latitude[currMeasNr] = np.nan

            navSolutions[0].longitude[currMeasNr] = np.nan

            navSolutions[0].height[currMeasNr] = np.nan

            navSolutions[0].E[currMeasNr] = np.nan

            navSolutions[0].N[currMeasNr] = np.nan

            navSolutions[0].U[currMeasNr] = np.nan

            channel[0].az[activeChnList, currMeasNr] = np.nan * np.ones(activeChnList.shape)

            channel[0].el[activeChnList, currMeasNr] = np.nan * np.ones(activeChnList.shape)

        # satellites are excluded do to elevation mask. Therefore rasing
        # satellites will be not included even if they will be above
        # elevation mask at some point. This would be a good place to
        # update positions of the excluded satellites.
        # === Update the transmit time ("measurement time") ====================
        transmitTime = transmitTime + settings.navSolPeriod / 1000
    return navSolutions, eph


if __name__ == '__main__':
    trackResults = np.load('trackingResults_python.npy')
    postNavigation(trackResults, initSettings.Settings())
