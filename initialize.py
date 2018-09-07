# ./initSettings.m

# Functions initializes and saves settings. Settings can be edited inside of
# the function, updated from the command line or updated using a dedicated
# GUI - "setSettings".

# All settings are described inside function code.

# settings = initSettings()

#   Inputs: none

#   Outputs:
#       settings     - Receiver settings (a structure).
import datetime

import numpy as np


class Result(object):
    def __init__(self, settings):
        self._settings = settings
        self._results = None
        self._channels = None

    @property
    def settings(self):
        return self._settings

    @property
    def channels(self):
        assert isinstance(self._channels, np.recarray)
        return self._channels

    @property
    def results(self):
        assert isinstance(self._results, np.recarray)
        return self._results

    @results.setter
    def results(self, records):
        assert isinstance(records, np.recarray)
        self._results = records

    def plot(self):
        pass


class TruePosition(object):
    def __init__(self):
        self._E = None
        self._N = None
        self._U = None

    @property
    def E(self):
        return self._E

    @E.setter
    def E(self, e):
        self._E = e

    @property
    def N(self):
        return self._N

    @N.setter
    def N(self, n):
        self._N = n

    @property
    def U(self):
        return self._U

    @U.setter
    def U(self, u):
        self._U = u


class Settings(object):
    def __init__(self):
        # Processing settings ====================================================
        # Number of milliseconds to be processed used 36000 + any transients (see
        # below - in Nav parameters) to ensure nav subframes are provided
        self.msToProcess = 37000.0

        # Number of channels to be used for signal processing
        self.numberOfChannels = 8

        # Move the starting point of processing. Can be used to start the signal
        # processing at any point in the data record (e.g. for long records). fseek
        # function is used to move the file read point, therefore advance is byte
        # based only.
        self.skipNumberOfBytes = 0

        # Raw signal file name and other parameter ===============================
        # This is a "default" name of the data file (signal record) to be used in
        # the post-processing mode
        self.fileName = '/Users/yangsu/Downloads/GNSS_signal_records/GPSdata-DiscreteComponents-fs38_192-if9_55.bin'

        # Data type used to store one sample
        self.dataType = 'int8'

        # Intermediate, sampling and code frequencies
        self.IF = 9548000.0

        self.samplingFreq = 38192000.0

        self.codeFreqBasis = 1023000.0

        # Define number of chips in a code period
        self.codeLength = 1023

        # Acquisition settings ===================================================
        # Skips acquisition in the script postProcessing.m if set to 1
        self.skipAcquisition = False

        # List of satellites to look for. Some satellites can be excluded to speed
        # up acquisition
        self.acqSatelliteList = range(1, 33)

        # Band around IF to search for satellite signal. Depends on max Doppler
        self.acqSearchBand = 14.0

        # Threshold for the signal presence decision rule
        self.acqThreshold = 2.5

        # Tracking loops settings ================================================
        # Code tracking loop parameters
        self.dllDampingRatio = 0.7

        self.dllNoiseBandwidth = 2.0

        self.dllCorrelatorSpacing = 0.5

        # Carrier tracking loop parameters
        self.pllDampingRatio = 0.7

        self.pllNoiseBandwidth = 25.0

        # Navigation solution settings ===========================================

        # Period for calculating pseudoranges and position
        self.navSolPeriod = 500.0

        # Elevation mask to exclude signals from satellites at low elevation
        self.elevationMask = 10.0

        # Enable/dissable use of tropospheric correction
        self.useTropCorr = True

        # 1 - On

        # True position of the antenna in UTM system (if known). Otherwise enter
        # all NaN's and mean position will be used as a reference .
        self.truePosition = TruePosition()
        #         self.truePosition.E = np.nan

        #         self.truePosition.N = np.nan

        #         self.truePosition.U = np.nan

        # Plot settings ==========================================================
        # Enable/disable plotting of the tracking results for each channel
        self.plotTracking = True

        # 1 - On

        # Constants ==============================================================

        self._c = 299792458.0

        self._startOffset = 68.802

    @property
    def c(self):
        return self._c

    @property
    def startOffset(self):
        return self._startOffset

    @property
    def samplesPerCode(self):
        return np.long(np.round(self.samplingFreq / (self.codeFreqBasis / self.codeLength)))

    # makeCaTable.m
    def makeCaTable(self):
        # Function generates CA codes for all 32 satellites based on the settings
        # provided in the structure "settings". The codes are digitized at the
        # sampling frequency specified in the settings structure.
        # One row in the "caCodesTable" is one C/A code. The row number is the PRN
        # number of the C/A code.

        # caCodesTable = makeCaTable(settings)

        #   Inputs:
        #       settings        - receiver settings
        #   Outputs:
        #       caCodesTable    - an array of arrays (matrix) containing C/A codes
        #                       for all satellite PRN-s

        # --- Find number of samples per spreading code ----------------------------
        samplesPerCode = self.samplesPerCode

        # --- Prepare the output matrix to speed up function -----------------------
        caCodesTable = np.zeros((32, samplesPerCode))

        # --- Find time constants --------------------------------------------------
        ts = 1.0 / self.samplingFreq

        tc = 1.0 / self.codeFreqBasis

        # === For all satellite PRN-s ...
        for PRN in range(32):
            # --- Generate CA code for given PRN -----------------------------------
            caCode = self.generateCAcode(PRN)

            # --- Make index array to read C/A code values -------------------------
            # The length of the index array depends on the sampling frequency -
            # number of samples per millisecond (because one C/A code period is one
            # millisecond).
            codeValueIndex = np.ceil(ts * np.arange(1, samplesPerCode + 1) / tc) - 1
            codeValueIndex = np.longlong(codeValueIndex)

            codeValueIndex[-1] = 1022

            # The "upsampled" code is made by selecting values form the CA code
            # chip array (caCode) for the time instances of each sample.
            caCodesTable[PRN] = caCode[codeValueIndex]
        return caCodesTable

    # generateCAcode.m
    def generateCAcode(self, prn):
        # generateCAcode.m generates one of the 32 GPS satellite C/A codes.

        # CAcode = generateCAcode(PRN)

        #   Inputs:
        #       PRN         - PRN number of the sequence.

        #   Outputs:
        #       CAcode      - a vector containing the desired C/A code sequence
        #                   (chips).

        # --- Make the code shift array. The shift depends on the PRN number -------
        # The g2s vector holds the appropriate shift of the g2 code to generate
        # the C/A code (ex. for SV#19 - use a G2 shift of g2s(19) = 471)

        assert prn in range(0, 32)
        g2s = [5, 6, 7, 8, 17, 18, 139, 140, 141, 251,
               252, 254, 255, 256, 257, 258, 469, 470, 471, 472,
               473, 474, 509, 512, 513, 514, 515, 516, 859, 860,
               861, 862,
               145, 175, 52, 21, 237, 235, 886, 657, 634, 762, 355, 1012, 176, 603, 130, 359, 595, 68, 386]

        # --- Pick right shift for the given PRN number ----------------------------
        g2shift = g2s[prn]

        # --- Generate G1 code -----------------------------------------------------

        # --- Initialize g1 output to speed up the function ---
        g1 = np.zeros(1023)

        # --- Load shift register ---
        reg = -1 * np.ones(10)

        # --- Generate all G1 signal chips based on the G1 feedback polynomial -----
        for i in range(1023):
            g1[i] = reg[-1]

            saveBit = reg[2] * reg[9]

            reg[1:] = reg[:-1]

            reg[0] = saveBit

        # --- Generate G2 code -----------------------------------------------------

        # --- Initialize g2 output to speed up the function ---
        g2 = np.zeros(1023)

        # --- Load shift register ---
        reg = -1 * np.ones(10)

        # --- Generate all G2 signal chips based on the G2 feedback polynomial -----
        for i in range(1023):
            g2[i] = reg[-1]

            saveBit = reg[1] * reg[2] * reg[5] * reg[7] * reg[8] * reg[9]

            reg[1:] = reg[:-1]

            reg[0] = saveBit

        # --- Shift G2 code --------------------------------------------------------
        # The idea: g2 = concatenate[ g2_right_part, g2_left_part ];
        g2 = np.r_[g2[1023 - g2shift:], g2[:1023 - g2shift]]

        # --- Form single sample C/A code by multiplying G1 and G2 -----------------
        CAcode = -g1 * g2
        return CAcode

    @staticmethod
    # calcLoopCoef.m
    def calcLoopCoef(LBW, zeta, k):
        # Function finds loop coefficients. The coefficients are used then in PLL-s
        # and DLL-s.

        # [tau1, tau2] = calcLoopCoef(LBW, zeta, k)

        #   Inputs:
        #       LBW           - Loop noise bandwidth
        #       zeta          - Damping ratio
        #       k             - Loop gain

        #   Outputs:
        #       tau1, tau2   - Loop filter coefficients

        # Solve natural frequency
        Wn = LBW * 8.0 * zeta / (4.0 * zeta ** 2 + 1)

        # solve for t1 & t2
        tau1 = k / (Wn * Wn)

        tau2 = 2.0 * zeta / Wn

        return tau1, tau2

    def probeData(self, fileNameStr=None):

        import matplotlib.pyplot as plt
        from scipy.signal import welch
        from scipy.signal.windows.windows import hamming

        # Function plots raw data information: time domain plot, a frequency domain
        # plot and a histogram.

        # The function can be called in two ways:
        #   probeData(settings)
        # or
        #   probeData(fileName, settings)

        #   Inputs:
        #       fileName        - name of the data file. File name is read from
        #                       settings if parameter fileName is not provided.

        #       settings        - receiver settings. Type of data file, sampling
        #                       frequency and the default filename are specified
        #                       here.

        # Check the number of arguments ==========================================
        if fileNameStr is None:
            fileNameStr = self.fileName
        if not isinstance(fileNameStr, str):
            raise TypeError('File name must be a string')
        settings = self
        # Generate plot of raw data ==============================================

        try:
            with open(fileNameStr, 'rb') as fid:
                # Move the starting point of processing. Can be used to start the
                # signal processing at any point in the data record (e.g. for long
                # records).
                fid.seek(settings.skipNumberOfBytes, 0)
                samplesPerCode = settings.samplesPerCode

                try:
                    data = np.fromfile(fid,
                                       settings.dataType,
                                       10 * samplesPerCode)

                except IOError:
                    # The file is too short
                    print 'Could not read enough data from the data file.'
                # --- Initialization ---------------------------------------------------
                plt.figure(100)
                plt.clf()
                timeScale = np.arange(0, 0.005, 1 / settings.samplingFreq)

                plt.subplot(2, 2, 1)
                plt.plot(1000 * timeScale[1:samplesPerCode / 50],
                         data[1:samplesPerCode / 50])
                plt.axis('tight')
                plt.grid()
                plt.title('Time domain plot')
                plt.xlabel('Time (ms)')
                plt.ylabel('Amplitude')
                plt.subplot(2, 2, 2)
                f, Pxx = welch(data - np.mean(data),
                               settings.samplingFreq / 1000000.0,
                               hamming(16384, False),
                               16384,
                               1024,
                               16384)
                plt.semilogy(f, Pxx)
                plt.axis('tight')
                plt.grid()
                plt.title('Frequency domain plot')
                plt.xlabel('Frequency (MHz)')
                plt.ylabel('Magnitude')
                plt.show()
                plt.subplot(2, 2, 3.5)
                plt.hist(data, np.arange(- 128, 128))
                dmax = np.max(np.abs(data)) + 1

                plt.axis('tight')
                adata = plt.axis()

                plt.axis([-dmax, dmax, adata[2], adata[3]])
                plt.grid('on')
                plt.title('Histogram')
                plt.xlabel('Bin')
                plt.ylabel('Number in bin')
            # === Error while opening the data file ================================
        except IOError as e:
            print 'Unable to read file "%s": %s' % (fileNameStr, e)

    # ./postProcessing.m

    # Script postProcessing.m processes the raw signal from the specified data
    # file (in settings) operating on blocks of 37 seconds of data.

    # First it runs acquisition code identifying the satellites in the file,
    # then the code and carrier for each of the satellites are tracked, storing
    # the 1m sec accumulations.  After processing all satellites in the 37 sec
    # data block, then postNavigation is called. It calculates pseudoranges
    # and attempts a position solutions. At the end plots are made for that
    # block of data.

    #                         THE SCRIPT "RECIPE"

    # The purpose of this script is to combine all parts of the software
    # receiver.

    # 1.1) Open the data file for the processing and seek to desired point.

    # 2.1) Acquire satellites

    # 3.1) Initialize channels (preRun.m).
    # 3.2) Pass the channel structure and the file identifier to the tracking
    # function. It will read and process the data. The tracking results are
    # stored in the trackResults structure. The results can be accessed this
    # way (the results are stored each millisecond):
    # trackResults(channelNumber).XXX(fromMillisecond : toMillisecond), where
    # XXX is a field name of the result (e.g. I_P, codePhase etc.)

    # 4) Pass tracking results to the navigation solution function. It will
    # decode navigation messages, find satellite positions, measure
    # pseudoranges and find receiver position.

    # 5) Plot the results.

    def postProcessing(self, fileNameStr=None):
        # Initialization =========================================================
        import acquisition
        import postNavigation
        import tracking
        print 'Starting processing...'
        settings = self
        if not fileNameStr:
            fileNameStr = settings.fileName
        if not isinstance(fileNameStr, str):
            raise TypeError('File name must be a string')
        try:
            with open(fileNameStr, 'rb') as fid:

                # If success, then process the data
                # Move the starting point of processing. Can be used to start the
                # signal processing at any point in the data record (e.g. good for long
                # records or for signal processing in blocks).
                fid.seek(settings.skipNumberOfBytes, 0)
                # Acquisition ============================================================
                # Do acquisition if it is not disabled in settings or if the variable
                # acqResults does not exist.
                if not settings.skipAcquisition:  # or 'acqResults' not in globals():
                    # Find number of samples per spreading code
                    samplesPerCode = settings.samplesPerCode

                    # frequency estimation
                    data = np.fromfile(fid, settings.dataType, 11 * samplesPerCode)

                    print '   Acquiring satellites...'
                    acqResults = acquisition.AcquisitionResult(settings)
                    acqResults.acquire(data)
                    acqResults.plot()
                # Initialize channels and prepare for the run ============================
                # Start further processing only if a GNSS signal was acquired (the
                # field FREQUENCY will be set to 0 for all not acquired signals)
                if np.any(acqResults.carrFreq):
                    acqResults.preRun()
                    acqResults.showChannelStatus()
                else:
                    # No satellites to track, exit
                    print 'No GNSS signals detected, signal processing finished.'
                    trackResults = None

                # Track the signal =======================================================
                startTime = datetime.datetime.now()

                print '   Tracking started at %s' % startTime.strftime('%X')
                trackResults = tracking.TrackingResult(acqResults)
                try:
                    trackResults.results = np.load('trackingResults_python.npy')
                except IOError:
                    trackResults.track(fid)
                    np.save('trackingResults_python', trackResults.results)

                print '   Tracking is over (elapsed time %s s)' % (datetime.datetime.now() - startTime).total_seconds()
                # Auto save the acquisition & tracking results to save time.
                print '   Saving Acquisition & Tracking results to storage'
                # Calculate navigation solutions =========================================
                print '   Calculating navigation solutions...'
                navResults = postNavigation.NavigationResult(trackResults)
                navResults.postNavigate()

                print '   Processing is complete for this data block'
                # Plot all results ===================================================
                print '   Plotting results...'
                # TODO turn off tracking plots for now
                if not settings.plotTracking:
                    trackResults.plot()
                navResults.plot()
                print 'Post processing of the signal is over.'
        except IOError as e:
            # Error while opening the data file.
            print 'Unable to read file "%s": %s.' % (settings.fileName, e)
