import numpy as np

import ephemeris
from geoFunctions import satpos, leastSquarePos, cart2geo, findUtmZone, cart2utm
from initialize import Result


class NavigationResult(Result):
    def __init__(self, trackResult):
        self._results = trackResult.results
        self._channels = trackResult.channels
        self._settings = trackResult.settings
        self._solutions = None
        self._eph = None

    @property
    def solutions(self):
        assert isinstance(self._solutions, np.recarray)
        return self._solutions

    @property
    def ephemeris(self):
        assert isinstance(self._solutions, np.recarray)
        return self._eph

    # ./calculatePseudoranges.m
    def calculatePseudoranges(self, msOfTheSignal, channelList):
        trackResults = self._results
        settings = self._settings
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
        samplesPerCode = settings.samplesPerCode

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
    def postNavigate(self):
        trackResults = self._results
        settings = self._settings
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

        # Check is there enough data to obtain any navigation solution ===========
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
            navSolutions = None
            self._solutions = navSolutions
            eph = None
            self._eph = eph
            return

        # Find preamble start positions ==========================================

        subFrameStart, activeChnList = self.findPreambles()

        # Decode ephemerides =====================================================
        field_str = 'weekNumber,accuracy,health,T_GD,IODC,t_oc,a_f2,a_f1,a_f0,'
        field_str += 'IODE_sf2,C_rs,deltan,M_0,C_uc,e,C_us,sqrtA,t_oe,'
        field_str += 'C_ic,omega_0,C_is,i_0,C_rc,omega,omegaDot,IODE_sf3,iDot'
        eph = np.recarray((32,), formats=['O'] * 27, names=field_str)
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

        # Check if the number of satellites is still above 3 =====================
        if activeChnList.size == 0 or activeChnList.size < 4:
            # Show error message and exit
            print 'Too few satellites with ephemeris data for position calculations. Exiting!'
            navSolutions = None
            self._solutions = navSolutions
            eph = None
            self._eph = eph
            return

        # Initialization =========================================================

        # Set the satellite elevations array to INF to include all satellites for
        # the first calculation of receiver position. There is no reference point
        # to find the elevation angle as there is no receiver position estimate at
        # this point.
        satElev = np.Inf * np.ones(settings.numberOfChannels)

        # Save the active channel list. The list contains satellites that are
        # tracked and have the required ephemeris data. In the next step the list
        # will depend on each satellite's elevation angle, which will change over
        # time.
        readyChnList = activeChnList.copy()

        transmitTime = TOW

        ###########################################################################
        #   Do the satellite and receiver position calculations                  #
        ###########################################################################
        # Initialization of current measurement ==================================
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
            # Exclude satellites, that are below elevation mask
            activeChnList = np.intersect1d((satElev >= settings.elevationMask).nonzero()[0], readyChnList)

            channel[0].PRN[activeChnList, currMeasNr] = trackResults[activeChnList].PRN

            # do to elevation mask will not "jump" to position (0,0) in the sky
            # plot.
            # channel[0].el[:, currMeasNr] = np.nan * np.ones(settings.numberOfChannels)

            # channel[0].az[:, currMeasNr] = np.nan * np.ones(settings.numberOfChannels)

            # Find pseudoranges ======================================================
            channel[0].rawP[:, currMeasNr] = self.calculatePseudoranges(
                subFrameStart + settings.navSolPeriod * currMeasNr,
                activeChnList)

            # Find satellites positions and clocks corrections =======================
            satPositions, satClkCorr = satpos(transmitTime, trackResults[activeChnList].PRN, eph, settings)

            # Find receiver position =================================================
            # 3D receiver position can be found only if signals from more than 3
            # satellites are available
            if activeChnList.size > 3:
                # === Calculate receiver position ==================================
                (xyzdt,
                 channel[0].el[activeChnList, currMeasNr],
                 channel[0].az[activeChnList, currMeasNr],
                 navSolutions[0].DOP[:, currMeasNr]) = leastSquarePos(satPositions,
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

                # Coordinate conversion ==================================================
                # === Convert to geodetic coordinates ==============================
                (navSolutions[0].latitude[currMeasNr],
                 navSolutions[0].longitude[currMeasNr],
                 navSolutions[0].height[currMeasNr]) = cart2geo(navSolutions[0].X[currMeasNr],
                                                                navSolutions[0].Y[currMeasNr],
                                                                navSolutions[0].Z[currMeasNr],
                                                                4)

                navSolutions[0].utmZone = findUtmZone(navSolutions[0].latitude[currMeasNr],
                                                      navSolutions[0].longitude[currMeasNr])

                (navSolutions[0].E[currMeasNr],
                 navSolutions[0].N[currMeasNr],
                 navSolutions[0].U[currMeasNr]) = cart2utm(xyzdt[0], xyzdt[1], xyzdt[2],
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

            # satellites are excluded do to elevation mask. Therefore raising
            # satellites will be not included even if they will be above
            # elevation mask at some point. This would be a good place to
            # update positions of the excluded satellites.
            # === Update the transmit time ("measurement time") ====================
            transmitTime += settings.navSolPeriod / 1000

        self._solutions = navSolutions
        self._eph = eph
        return

    def plot(self):
        settings = self._settings
        navSolutions = self._solutions
        assert isinstance(navSolutions, np.recarray)

        import matplotlib as mpl
        import matplotlib.gridspec as gs
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import axes3d

        import initialize

        # %% configure matplotlib
        mpl.rcdefaults()
        # mpl.rcParams['font.sans-serif']
        # mpl.rcParams['font.family'] = 'serif'
        mpl.rc('savefig', bbox='tight', transparent=False, format='png')
        mpl.rc('axes', grid=True, linewidth=1.5, axisbelow=True)
        mpl.rc('lines', linewidth=1.5, solid_joinstyle='bevel')
        mpl.rc('figure', figsize=[8, 6], autolayout=False, dpi=120)
        mpl.rc('text', usetex=True)
        mpl.rc('font', family='serif', serif='Computer Modern Roman', size=10)
        mpl.rc('mathtext', fontset='cm')

        # mpl.rc('font', size=16)
        # mpl.rc('text.latex', preamble=r'\usepackage{cmbright}')

        # ./plotNavigation.m

        # Functions plots variations of coordinates over time and a 3D position
        # plot. It plots receiver coordinates in UTM system or coordinate offsets if
        # the true UTM receiver coordinates are provided.

        # plotNavigation(navSolutions, settings)

        #   Inputs:
        #       navSolutions    - Results from navigation solution function. It
        #                       contains measured pseudoranges and receiver
        #                       coordinates.
        #       settings        - Receiver settings. The true receiver coordinates
        #                       are contained in this structure.

        # Plot results in the necessary data exists ==============================
        if navSolutions is not None:
            refCoord = initialize.TruePosition()
            # If reference position is not provided, then set reference position
            # to the average postion
            if settings.truePosition.E is None or settings.truePosition.N is None or settings.truePosition.U is None:
                # === Compute mean values ==========================================
                # Remove NaN-s or the output of the function MEAN will be NaN.
                refCoord.E = np.nanmean(navSolutions[0].E)

                refCoord.N = np.nanmean(navSolutions[0].N)

                refCoord.U = np.nanmean(navSolutions[0].U)

                meanLongitude = np.nanmean(navSolutions[0].longitude)

                meanLatitude = np.nanmean(navSolutions[0].latitude)

                refPointLgText = 'Mean Position' + '\\newline Lat: %.5f $^\circ$' % meanLatitude + \
                                 '\\newline Lng: %.5f $^\circ$' % meanLongitude + \
                                 '\\newline Hgt: %+6.1f' % np.nanmean(navSolutions[0].height)

            else:
                refPointLgText = 'Reference Position'

                refCoord.E = settings.truePosition.E

                refCoord.N = settings.truePosition.N

                refCoord.U = settings.truePosition.U

            figureNumber = 300

            # figure windows, when many figures are closed and reopened. Figures
            # drawn or opened by the user, will not be "overwritten" by this
            # function if the auto numbering is not used.
            # === Select (or create) and clear the figure ==========================
            f = plt.figure(figureNumber)
            f.clf()
            f.set_label('Navigation solutions')
            spec = gs.GridSpec(2, 2)
            h11 = plt.subplot(spec[0:2])

            # the axes3d module is needed for the following line
            dummy = axes3d.Axes3D
            h31 = plt.subplot(spec[2], projection='3d')

            h32 = plt.subplot(spec[3], projection='polar')

            # Plot all figures =======================================================
            # --- Coordinate differences in UTM system -----------------------------
            h11.plot(navSolutions[0].E - refCoord.E, '-',
                     navSolutions[0].N - refCoord.N, '-',
                     navSolutions[0].U - refCoord.U, '-')
            h11.legend(['E', 'N', 'U'])
            h11.set(title='Coordinates variations in UTM system',
                    xlabel='Measurement period: %i ms' % settings.navSolPeriod,
                    ylabel='Variations (m)')
            h11.grid()
            h11.axis('tight')
            h31.plot((navSolutions[0].E - refCoord.E).T,
                     (navSolutions[0].N - refCoord.N).T,
                     (navSolutions[0].U - refCoord.U).T, '+')
            h31.hold(True)
            h31.plot([0], [0], [0], 'r+', lw=1.5, ms=10)
            h31.hold(False)
            # h31.viewLim(0,90)
            h31.axis('equal')
            h31.grid(which='minor')
            h31.legend(['Measurements', refPointLgText])
            h31.set(title='Positions in UTM system (3D plot)',
                    xlabel='East (m)',
                    ylabel='North (m)',
                    zlabel='Upping (m)')
            h32.plot(np.deg2rad(navSolutions[0].channel[0].az.T),
                     90 - navSolutions[0].channel[0].el.T)
            [h32.text(x, y, s) for x, y, s in zip(np.deg2rad(navSolutions[0].channel[0].az[:, 0]),
                                                  90 - navSolutions[0].channel[0].el[:, 0],
                                                  navSolutions[0].channel[0].PRN[:, 0])]
            h32.set_theta_direction(-1)
            h32.set_theta_zero_location('N')
            h32.set_xlim([0, 2 * np.pi])
            h32.set_xticks(np.linspace(0, 2 * np.pi, 12, endpoint=False))
            h32.set_rlabel_position(0)
            h32.set_ylim([0, 90])
            h32.set_yticks([0, 15, 30, 45, 60, 75])
            h32.set_yticklabels([90, 75, 60, 45, 30, 15])
            h32.set_title('Sky plot (mean PDOP: %f )' % np.mean(navSolutions[0].DOP[1, :]))
            f.show()
        else:
            print 'plotNavigation: No navigation data to plot.'

    @staticmethod
    # navPartyChk.m
    def navPartyChk(ndat):
        # This function is called to compute and status the parity bits on GPS word.
        # Based on the flowchart in Figure 2-10 in the 2nd Edition of the GPS-SPS
        # Signal Spec.

        # status = navPartyChk(ndat)

        #   Inputs:
        #       ndat        - an array (1x32) of 32 bits represent a GPS navigation
        #                   word which is 30 bits plus two previous bits used in
        #                   the parity calculation (-2 -1 0 1 2 ... 28 29)

        #   Outputs:
        #       status      - the test value which equals EITHER +1 or -1 if parity
        #                   PASSED or 0 if parity fails.  The +1 means bits #1-24
        #                   of the current word have the correct polarity, while -1
        #                   means the bits #1-24 of the current word must be
        #                   inverted.

        # In order to accomplish the exclusive or operation using multiplication
        # this program represents a '0' with a '-1' and a '1' with a '1' so that
        # the exclusive or table holds true for common data operations

        #	a	b	xor 			a	b	product
        #  --------------          -----------------
        #	0	0	 1			   -1  -1	   1
        #	0	1	 0			   -1   1	  -1
        #	1	0	 0			    1  -1	  -1
        #	1	1	 1			    1   1	   1

        # --- Check if the data bits must be inverted ------------------------------
        if ndat[1] != 1:
            ndat[2:26] *= (-1)

        # --- Calculate 6 parity bits ----------------------------------------------
        # The elements of the ndat array correspond to the bits showed in the table
        # 20-XIV (ICD-200C document) in the following way:
        # The first element in the ndat is the D29* bit and the second - D30*.
        # The elements 3 - 26 are bits d1-d24 in the table.
        # The elements 27 - 32 in the ndat array are the received bits D25-D30.
        # The array "parity" contains the computed D25-D30 (parity) bits.
        parity = np.zeros(6)
        parity[0] = ndat[0] * ndat[2] * ndat[3] * ndat[4] * ndat[6] * \
                    ndat[7] * ndat[11] * ndat[12] * ndat[13] * ndat[14] * \
                    ndat[15] * ndat[18] * ndat[19] * ndat[21] * ndat[24]

        parity[1] = ndat[1] * ndat[3] * ndat[4] * ndat[5] * ndat[7] * \
                    ndat[8] * ndat[12] * ndat[13] * ndat[14] * ndat[15] * \
                    ndat[16] * ndat[19] * ndat[20] * ndat[22] * ndat[25]

        parity[2] = ndat[0] * ndat[2] * ndat[4] * ndat[5] * ndat[6] * \
                    ndat[8] * ndat[9] * ndat[13] * ndat[14] * ndat[15] * \
                    ndat[16] * ndat[17] * ndat[20] * ndat[21] * ndat[23]

        parity[3] = ndat[1] * ndat[3] * ndat[5] * ndat[6] * ndat[7] * \
                    ndat[9] * ndat[10] * ndat[14] * ndat[15] * ndat[16] * \
                    ndat[17] * ndat[18] * ndat[21] * ndat[22] * ndat[24]

        parity[4] = ndat[1] * ndat[2] * ndat[4] * ndat[6] * ndat[7] * \
                    ndat[8] * ndat[10] * ndat[11] * ndat[15] * ndat[16] * \
                    ndat[17] * ndat[18] * ndat[19] * ndat[22] * ndat[23] * \
                    ndat[25]

        parity[5] = ndat[0] * ndat[4] * ndat[6] * ndat[7] * ndat[9] * \
                    ndat[10] * ndat[11] * ndat[12] * ndat[14] * ndat[16] * \
                    ndat[20] * ndat[23] * ndat[24] * ndat[25]

        # --- Compare if the received parity is equal the calculated parity --------
        if (parity == ndat[26:]).sum() == 6:
            # Parity is OK. Function output is -1 or 1 depending if the data bits
            # must be inverted or not. The "ndat[2]" is D30* bit - the last  bit of
            # previous subframe.
            status = -1 * ndat[1]

        else:
            # Parity failure
            status = 0

        return status

    # ./findPreambles.m
    def findPreambles(self):
        assert isinstance(self._results, np.recarray)
        trackResults = self._results
        settings = self._settings
        # findPreambles finds the first preamble occurrence in the bit stream of
        # each channel. The preamble is verified by check of the spacing between
        # preambles (6sec) and parity checking of the first two words in a
        # subframe. At the same time function returns list of channels, that are in
        # tracking state and with valid preambles in the nav data stream.

        # [firstSubFrame, activeChnList] = findPreambles(trackResults, settings)

        #   Inputs:
        #       trackResults    - output from the tracking function
        #       settings        - Receiver settings.

        #   Outputs:
        #       firstSubframe   - the array contains positions of the first
        #                       preamble in each channel. The position is ms count
        #                       since start of tracking. Corresponding value will
        #                       be set to 0 if no valid preambles were detected in
        #                       the channel.
        #       activeChnList   - list of channels containing valid preambles

        # Preamble search can be delayed to a later point in the tracking results
        # to avoid noise due to tracking loop transients
        searchStartOffset = 0

        # --- Initialize the firstSubFrame array -----------------------------------
        firstSubFrame = np.zeros(settings.numberOfChannels, dtype=int)

        # --- Generate the preamble pattern ----------------------------------------
        preamble_bits = np.r_[1, - 1, - 1, - 1, 1, - 1, 1, 1]

        # "Upsample" the preamble - make 20 vales per one bit. The preamble must be
        # found with precision of a sample.
        preamble_ms = np.kron(preamble_bits, np.ones(20))

        # --- Make a list of channels excluding not tracking channels --------------
        activeChnList = (trackResults.status != '-').nonzero()[0]

        # === For all tracking channels ...
        for channelNr in range(len(activeChnList)):
            # Correlate tracking output with preamble ================================
            # Read output from tracking. It contains the navigation bits. The start
            # of record is skipped here to avoid tracking loop transients.
            bits = trackResults[channelNr].I_P[searchStartOffset:].copy()

            bits[bits > 0] = 1

            bits[bits <= 0] = - 1

            # have to zero pad the preamble so that they are the same length
            tlmXcorrResult = np.correlate(bits,
                                          np.pad(preamble_ms, (0, bits.size - preamble_ms.size), 'constant'),
                                          mode='full')

            # Find all starting points off all preamble like patterns ================
            # clear('index')
            # clear('index2')
            xcorrLength = (len(tlmXcorrResult) + 1) / 2

            index = (np.abs(tlmXcorrResult[xcorrLength - 1:xcorrLength * 2]) > 153).nonzero()[0] + searchStartOffset

            # Analyze detected preamble like patterns ================================
            for i in range(len(index)):
                # --- Find distances in time between this occurrence and the rest of
                # preambles like patterns. If the distance is 6000 milliseconds (one
                # subframe), the do further verifications by validating the parities
                # of two GPS words
                index2 = index - index[i]

                if (index2 == 6000).any():
                    # === Re-read bit vales for preamble verification ==============
                    # Preamble occurrence is verified by checking the parity of
                    # the first two words in the subframe. Now it is assumed that
                    # bit boundaries a known. Therefore the bit values over 20ms are
                    # combined to increase receiver performance for noisy signals.
                    # in Total 62 bits mast be read :
                    # 2 bits from previous subframe are needed for parity checking;
                    # 60 bits for the first two 30bit words (TLM and HOW words).
                    # The index is pointing at the start of TLM word.
                    bits = trackResults[channelNr].I_P[index[i] - 40:index[i] + 20 * 60].copy()

                    bits = bits.reshape(20, -1, order='F')

                    bits = bits.sum(0)

                    bits[bits > 0] = 1

                    bits[bits <= 0] = - 1

                    if self.navPartyChk(bits[:32]) != 0 and self.navPartyChk(bits[30:62]) != 0:
                        # Parity was OK. Record the preamble start position. Skip
                        # the rest of preamble pattern checking for this channel
                        # and process next channel.
                        firstSubFrame[channelNr] = index[i]

                        break
            # Exclude channel from the active channel list if no valid preamble was
            # detected
            if firstSubFrame[channelNr] == 0:
                # Exclude channel from further processing. It does not contain any
                # valid preamble and therefore nothing more can be done for it.
                activeChnList = np.setdiff1d(activeChnList, channelNr)

                print 'Could not find valid preambles in channel %2d !' % channelNr
        return firstSubFrame, activeChnList


if __name__ == '__main__':
    pass
