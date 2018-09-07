import numpy as np

from initialize import Result


class AcquisitionResult(Result):
    def __init__(self, settings):
        self._settings = settings
        self._results = None
        self._channels = None

    @property
    def peakMetric(self):
        assert isinstance(self._results, np.recarray)
        return self._results.peakMetric

    @property
    def carrFreq(self):
        assert isinstance(self._results, np.recarray)
        return self._results.carrFreq

    @property
    def codePhase(self):
        assert isinstance(self._results, np.recarray)
        return self._results.codePhase

    def acquire(self, longSignal):
        # ./acquisition.m
        # Function performs cold start acquisition on the collected "data". It
        # searches for GPS signals of all satellites, which are listed in field
        # "acqSatelliteList" in the settings structure. Function saves code phase
        # and frequency of the detected signals in the "acqResults" structure.

        # acqResults = acquisition(longSignal, settings)

        #   Inputs:
        #       longSignal    - 11 ms of raw signal from the front-end
        #       settings      - Receiver settings. Provides information about
        #                       sampling and intermediate frequencies and other
        #                       parameters including the list of the satellites to
        #                       be acquired.
        #   Outputs:
        #       acqResults    - Function saves code phases and frequencies of the
        #                       detected signals in the "acqResults" structure. The
        #                       field "carrFreq" is set to 0 if the signal is not
        #                       detected for the given PRN number.

        # Initialization =========================================================
        settings = self._settings

        # Find number of samples per spreading code
        samplesPerCode = settings.samplesPerCode

        # Create two 1m sec vectors of data to correlate with and one with zero DC
        signal1 = longSignal[0:samplesPerCode]

        signal2 = longSignal[samplesPerCode:2 * samplesPerCode]

        signal0DC = longSignal - longSignal.mean()

        # Find sampling period
        ts = 1.0 / settings.samplingFreq

        # Find phase points of the local carrier wave
        phasePoints = np.arange(samplesPerCode) * 2 * np.pi * ts

        # Number of the frequency bins for the given acquisition band (500Hz steps)
        numberOfFrqBins = np.int(np.round(settings.acqSearchBand * 2) + 1)

        # Generate all C/A codes and sample them according to the sampling freq.
        caCodesTable = settings.makeCaTable()

        # --- Initialize arrays to speed up the code -------------------------------
        # Search results of all frequency bins and code shifts (for one satellite)
        results = np.zeros((numberOfFrqBins, samplesPerCode))

        # Carrier frequencies of the frequency bins
        frqBins = np.zeros(numberOfFrqBins)

        # --- Initialize acqResults ------------------------------------------------
        # Carrier frequencies of detected signals
        carrFreq = np.zeros(32)

        # C/A code phases of detected signals
        codePhase_ = np.zeros(32)

        # Correlation peak ratios of the detected signals
        peakMetric = np.zeros(32)

        print '('
        # Perform search for all listed PRN numbers ...
        for PRN in range(len(settings.acqSatelliteList)):
            # Correlate signals ======================================================
            # --- Perform DFT of C/A code ------------------------------------------
            caCodeFreqDom = np.fft.fft(caCodesTable[PRN, :]).conj()

            for frqBinIndex in range(numberOfFrqBins):
                # --- Generate carrier wave frequency grid (0.5kHz step) -----------
                frqBins[frqBinIndex] = settings.IF - \
                                       settings.acqSearchBand / 2 * 1000 + \
                                       500.0 * frqBinIndex

                sinCarr = np.sin(frqBins[frqBinIndex] * phasePoints)

                cosCarr = np.cos(frqBins[frqBinIndex] * phasePoints)

                I1 = sinCarr * signal1

                Q1 = cosCarr * signal1

                I2 = sinCarr * signal2

                Q2 = cosCarr * signal2

                IQfreqDom1 = np.fft.fft(I1 + 1j * Q1)

                IQfreqDom2 = np.fft.fft(I2 + 1j * Q2)

                # domain)
                convCodeIQ1 = IQfreqDom1 * caCodeFreqDom

                convCodeIQ2 = IQfreqDom2 * caCodeFreqDom

                acqRes1 = abs(np.fft.ifft(convCodeIQ1)) ** 2

                acqRes2 = abs(np.fft.ifft(convCodeIQ2)) ** 2

                # "blend" 1st and 2nd msec but will correct data bit issues
                if acqRes1.max() > acqRes2.max():
                    results[frqBinIndex, :] = acqRes1

                else:
                    results[frqBinIndex, :] = acqRes2

            # Look for correlation peaks in the results ==============================
            # Find the highest peak and compare it to the second highest peak
            # The second peak is chosen not closer than 1 chip to the highest peak
            # --- Find the correlation peak and the carrier frequency --------------
            peakSize = results.max(1).max()
            frequencyBinIndex = results.max(1).argmax()

            peakSize = results.max(0).max()
            codePhase = results.max(0).argmax()

            samplesPerCodeChip = long(round(settings.samplingFreq / settings.codeFreqBasis))

            excludeRangeIndex1 = codePhase - samplesPerCodeChip

            excludeRangeIndex2 = codePhase + samplesPerCodeChip

            # boundaries
            if excludeRangeIndex1 <= 0:
                codePhaseRange = np.r_[excludeRangeIndex2:samplesPerCode + excludeRangeIndex1 + 1]

            elif excludeRangeIndex2 >= samplesPerCode - 1:
                codePhaseRange = np.r_[excludeRangeIndex2 - samplesPerCode:excludeRangeIndex1]

            else:
                codePhaseRange = np.r_[0:excludeRangeIndex1 + 1, excludeRangeIndex2:samplesPerCode]

            # --- Find the second highest correlation peak in the same freq. bin ---
            secondPeakSize = results[frequencyBinIndex, codePhaseRange].max()

            peakMetric[PRN] = peakSize / secondPeakSize

            if (peakSize / secondPeakSize) > settings.acqThreshold:
                # Fine resolution frequency search =======================================
                # --- Indicate PRN number of the detected signal -------------------
                print '%02d ' % (PRN + 1)
                caCode = settings.generateCAcode(PRN)

                codeValueIndex = np.floor(ts * np.arange(1, 10 * samplesPerCode + 1) / (1.0 / settings.codeFreqBasis))

                longCaCode = caCode[np.longlong(codeValueIndex % 1023)]

                # (Using detected C/A code phase)
                xCarrier = signal0DC[codePhase:codePhase + 10 * samplesPerCode] * longCaCode

                fftNumPts = 8 * 2 ** (np.ceil(np.log2(len(xCarrier))))

                # associated carrier frequency
                fftxc = np.abs(np.fft.fft(xCarrier, np.long(fftNumPts)))

                uniqFftPts = np.long(np.ceil((fftNumPts + 1) / 2.0))

                fftMax = fftxc[4:uniqFftPts - 5].max()
                fftMaxIndex = fftxc[4:uniqFftPts - 5].argmax()

                fftFreqBins = np.arange(uniqFftPts) * settings.samplingFreq / fftNumPts

                carrFreq[PRN] = fftFreqBins[fftMaxIndex]

                codePhase_[PRN] = codePhase

            else:
                # --- No signal with this PRN --------------------------------------
                print '. '

        # === Acquisition is over ==================================================
        print ')\n'
        acqResults = np.core.records.fromarrays([carrFreq, codePhase_, peakMetric],
                                                names='carrFreq,codePhase,peakMetric')
        self._results = acqResults
        return

    def plot(self):
        assert isinstance(self._results, np.recarray)
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        # from scipy.io.matlab import loadmat

        # %% configure matplotlib
        mpl.rcdefaults()
        # mpl.rcParams['font.sans-serif']
        # mpl.rcParams['font.family'] = 'serif'
        mpl.rc('savefig', bbox='tight', transparent=False, format='png')
        mpl.rc('axes', grid=True, linewidth=1.5, axisbelow=True)
        mpl.rc('lines', linewidth=1.5, solid_joinstyle='bevel')
        mpl.rc('figure', figsize=[8, 6], autolayout=False, dpi=120)
        mpl.rc('text', usetex=True)
        mpl.rc('font', family='serif', serif='Computer Modern Roman', size=16)
        mpl.rc('mathtext', fontset='cm')

        # mpl.rc('font', size=16)
        # mpl.rc('text.latex', preamble=r'\usepackage{cmbright}')

        # ./plotAcquisition.m
        # Functions plots bar plot of acquisition results (acquisition metrics). No
        # bars are shown for the satellites not included in the acquisition list (in
        # structure SETTINGS).

        # plotAcquisition(acqResults)

        #   Inputs:
        #       acqResults    - Acquisition results from function acquisition.

        # Plot all results =======================================================
        f, hAxes = plt.subplots()

        plt.bar(range(1, 33), self.peakMetric)
        plt.title('Acquisition results')
        plt.xlabel('PRN number (no bar - SV is not in the acquisition list)')
        plt.ylabel('Acquisition Metric ($1^{st}$ to $2^{nd}$ Correlation Peaks Ratio')
        oldAxis = plt.axis()

        plt.axis([0, 33, 0, oldAxis[-1]])
        plt.xticks(range(1, 33), size=12)
        # plt.minorticks_on()
        hAxes.xaxis.grid()
        # Mark acquired signals ==================================================

        acquiredSignals = self.peakMetric * (self.carrFreq > 0)

        plt.bar(range(1, 33), acquiredSignals, FaceColor=(0, 0.8, 0))
        plt.legend(['Not acquired signals', 'Acquired signals'])
        plt.show()

    # preRun.m
    def preRun(self):
        assert isinstance(self._results, np.recarray)
        # Function initializes tracking channels from acquisition data. The acquired
        # signals are sorted according to the signal strength. This function can be
        # modified to use other satellite selection algorithms or to introduce
        # acquired signal properties offsets for testing purposes.

        # [channel] = preRun(acqResults, settings)

        #   Inputs:
        #       acqResults  - results from acquisition.
        #       settings    - receiver settings

        #   Outputs:
        #       channel     - structure contains information for each channel (like
        #                   properties of the tracked signal, channel status etc.).

        settings = self._settings
        # Initialize all channels ================================================
        PRN = np.zeros(settings.numberOfChannels, dtype='int64')
        acquiredFreq = np.zeros(settings.numberOfChannels)
        codePhase = np.zeros(settings.numberOfChannels)
        status = ['-' for _ in range(settings.numberOfChannels)]

        # --- Copy initial data to all channels ------------------------------------

        # Copy acquisition results ===============================================

        # --- Sort peaks to find strongest signals, keep the peak index information
        PRNindexes = sorted(enumerate(self.peakMetric),
                            key=lambda x: x[-1], reverse=True)

        # --- Load information about each satellite --------------------------------
        # Maximum number of initialized channels is number of detected signals, but
        # not more as the number of channels specified in the settings.
        for ii in range(min(settings.numberOfChannels, sum(self.carrFreq > 0))):
            PRN[ii] = PRNindexes[ii][0] + 1

            acquiredFreq[ii] = self.carrFreq[PRNindexes[ii][0]]

            codePhase[ii] = self.codePhase[PRNindexes[ii][0]]

            status[ii] = 'T'

        channel = np.core.records.fromarrays([PRN, acquiredFreq, codePhase, status],
                                             names='PRN,acquiredFreq,codePhase,status')
        self._channels = channel
        return

    def showChannelStatus(self):
        # Prints the status of all channels in a table.

        # showChannelStatus(channel, settings)

        #   Inputs:
        #       channel     - data for each channel. It is used to initialize and
        #                   at the processing of the signal (tracking part).
        #       settings    - receiver settings

        channel = self._channels
        settings = self._settings
        assert isinstance(channel, np.recarray)
        print ('\n*=========*=====*===============*===========*=============*========*')
        print ('| Channel | PRN |   Frequency   |  Doppler  | Code Offset | Status |')
        print ('*=========*=====*===============*===========*=============*========*')
        for channelNr in range(settings.numberOfChannels):
            if channel[channelNr].status != '-':
                print '|      %2d | %3d |  %2.5e |   %5.0f   |    %6d   |     %1s  |' % (
                    channelNr,
                    channel[channelNr].PRN,
                    channel[channelNr].acquiredFreq,
                    channel[channelNr].acquiredFreq - settings.IF,
                    channel[channelNr].codePhase,
                    channel[channelNr].status)
            else:
                print '|      %2d | --- |  ------------ |   -----   |    ------   |   Off  |' % channelNr

        print '*=========*=====*===============*===========*=============*========*\n'


if __name__ == '__main__':
    pass
