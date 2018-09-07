import numpy as np

from initialize import Result


class TrackingResult(Result):
    def __init__(self, acqResult):
        self._results = None
        self._channels = acqResult.channels
        self._settings = acqResult.settings

    # ./tracking.m
    def track(self, fid):
        channel = self._channels
        settings = self._settings
        # Performs code and carrier tracking for all channels.

        # [trackResults, channel] = tracking(fid, channel, settings)

        #   Inputs:
        #       fid             - file identifier of the signal record.
        #       channel         - PRN, carrier frequencies and code phases of all
        #                       satellites to be tracked (prepared by preRum.m from
        #                       acquisition results).
        #       settings        - receiver settings.
        #   Outputs:
        #       trackResults    - tracking results (structure array). Contains
        #                       in-phase prompt outputs and absolute starting
        #                       positions of spreading codes, together with other
        #                       observation data from the tracking loops. All are
        #                       saved every millisecond.

        # Initialize tracking variables ==========================================

        codePeriods = settings.msToProcess

        # --- DLL variables --------------------------------------------------------
        # Define early-late offset (in chips)
        earlyLateSpc = settings.dllCorrelatorSpacing

        # Summation interval
        PDIcode = 0.001

        # Calculate filter coefficient values
        tau1code, tau2code = settings.calcLoopCoef(settings.dllNoiseBandwidth, settings.dllDampingRatio, 1.0)

        # --- PLL variables --------------------------------------------------------
        # Summation interval
        PDIcarr = 0.001

        # Calculate filter coefficient values
        tau1carr, tau2carr = settings.calcLoopCoef(settings.pllNoiseBandwidth, settings.pllDampingRatio, 0.25)

        #     hwb=waitbar(0,'Tracking...')

        # Initialize a temporary list of records
        rec = []
        # Start processing channels ==============================================
        for channelNr in range(settings.numberOfChannels):
            msToProcess = np.long(settings.msToProcess)
            # Initialize fields for record(structured) array of tracked results
            status = '-'

            # The absolute sample in the record of the C/A code start:
            absoluteSample = np.zeros(msToProcess)

            # Freq of the C/A code:
            codeFreq_ = np.Inf * np.ones(msToProcess)

            # Frequency of the tracked carrier wave:
            carrFreq_ = np.Inf * np.ones(msToProcess)

            # Outputs from the correlators (In-phase):
            I_P_ = np.zeros(msToProcess)

            I_E_ = np.zeros(msToProcess)

            I_L_ = np.zeros(msToProcess)

            # Outputs from the correlators (Quadrature-phase):
            Q_E_ = np.zeros(msToProcess)

            Q_P_ = np.zeros(msToProcess)

            Q_L_ = np.zeros(msToProcess)

            # Loop discriminators
            dllDiscr = np.Inf * np.ones(msToProcess)

            dllDiscrFilt = np.Inf * np.ones(msToProcess)

            pllDiscr = np.Inf * np.ones(msToProcess)

            pllDiscrFilt = np.Inf * np.ones(msToProcess)

            PRN = 0

            # Only process if PRN is non zero (acquisition was successful)
            if channel[channelNr].PRN != 0:
                # Save additional information - each channel's tracked PRN
                PRN = channel[channelNr].PRN

                # signal processing at any point in the data record (e.g. for long
                # records). In addition skip through that data file to start at the
                # appropriate sample (corresponding to code phase). Assumes sample
                # type is schar (or 1 byte per sample)
                fid.seek(settings.skipNumberOfBytes + channel[channelNr].codePhase, 0)
                # Here PRN is the actual satellite ID instead of the 0-based index
                caCode = settings.generateCAcode(channel[channelNr].PRN - 1)

                caCode = np.r_[caCode[-1], caCode, caCode[0]]

                # define initial code frequency basis of NCO
                codeFreq = settings.codeFreqBasis

                remCodePhase = 0.0

                carrFreq = channel[channelNr].acquiredFreq

                carrFreqBasis = channel[channelNr].acquiredFreq

                remCarrPhase = 0.0

                oldCodeNco = 0.0

                oldCodeError = 0.0

                oldCarrNco = 0.0

                oldCarrError = 0.0

                for loopCnt in range(np.long(codePeriods)):
                    # GUI update -------------------------------------------------------------
                    # The GUI is updated every 50ms. This way Matlab GUI is still
                    # responsive enough. At the same time Matlab is not occupied
                    # all the time with GUI task.
                    if loopCnt % 50 == 0:
                        try:
                            print 'Tracking: Ch %d' % (channelNr + 1) + ' of %d' % settings.numberOfChannels + \
                                  '; PRN#%02d' % channel[channelNr].PRN + \
                                  '; Completed %d' % loopCnt + ' of %d' % codePeriods + ' msec'
                        finally:
                            pass
                    # Read next block of data ------------------------------------------------
                    # Find the size of a "block" or code period in whole samples
                    # Update the phasestep based on code freq (variable) and
                    # sampling frequency (fixed)
                    codePhaseStep = codeFreq / settings.samplingFreq

                    blksize = np.ceil((settings.codeLength - remCodePhase) / codePhaseStep)
                    blksize = np.long(blksize)

                    # interaction
                    rawSignal = np.fromfile(fid, settings.dataType, blksize)
                    samplesRead = len(rawSignal)

                    # If did not read in enough samples, then could be out of
                    # data - better exit
                    if samplesRead != blksize:
                        print 'Not able to read the specified number of samples for tracking, exiting!'
                        fid.close()
                        trackResults = None
                        return trackResults
                    # Set up all the code phase tracking information -------------------------
                    # Define index into early code vector
                    tcode = np.linspace(remCodePhase - earlyLateSpc,
                                        blksize * codePhaseStep + remCodePhase - earlyLateSpc,
                                        blksize, endpoint=False)

                    tcode2 = np.ceil(tcode)

                    earlyCode = caCode[np.longlong(tcode2)]

                    tcode = np.linspace(remCodePhase + earlyLateSpc,
                                        blksize * codePhaseStep + remCodePhase + earlyLateSpc,
                                        blksize, endpoint=False)

                    tcode2 = np.ceil(tcode)

                    lateCode = caCode[np.longlong(tcode2)]

                    tcode = np.linspace(remCodePhase,
                                        blksize * codePhaseStep + remCodePhase,
                                        blksize, endpoint=False)

                    tcode2 = np.ceil(tcode)

                    promptCode = caCode[np.longlong(tcode2)]

                    remCodePhase = tcode[blksize - 1] + codePhaseStep - 1023.0

                    # Generate the carrier frequency to mix the signal to baseband -----------
                    time = np.arange(0, blksize + 1) / settings.samplingFreq

                    trigarg = carrFreq * 2.0 * np.pi * time + remCarrPhase

                    remCarrPhase = trigarg[blksize] % (2 * np.pi)

                    carrCos = np.cos(trigarg[0:blksize])

                    carrSin = np.sin(trigarg[0:blksize])

                    # Generate the six standard accumulated values ---------------------------
                    # First mix to baseband
                    qBasebandSignal = carrCos * rawSignal

                    iBasebandSignal = carrSin * rawSignal

                    I_E = (earlyCode * iBasebandSignal).sum()

                    Q_E = (earlyCode * qBasebandSignal).sum()

                    I_P = (promptCode * iBasebandSignal).sum()

                    Q_P = (promptCode * qBasebandSignal).sum()

                    I_L = (lateCode * iBasebandSignal).sum()

                    Q_L = (lateCode * qBasebandSignal).sum()

                    # Find PLL error and update carrier NCO ----------------------------------
                    # Implement carrier loop discriminator (phase detector)
                    carrError = np.arctan(Q_P / I_P) / 2.0 / np.pi

                    carrNco = oldCarrNco + \
                              tau2carr / tau1carr * (carrError - oldCarrError) + \
                              carrError * (PDIcarr / tau1carr)

                    oldCarrNco = carrNco

                    oldCarrError = carrError

                    carrFreq = carrFreqBasis + carrNco

                    carrFreq_[loopCnt] = carrFreq

                    # Find DLL error and update code NCO -------------------------------------
                    codeError = (np.sqrt(I_E * I_E + Q_E * Q_E) - np.sqrt(I_L * I_L + Q_L * Q_L)) / (
                            np.sqrt(I_E * I_E + Q_E * Q_E) + np.sqrt(I_L * I_L + Q_L * Q_L))

                    codeNco = oldCodeNco + \
                              tau2code / tau1code * (codeError - oldCodeError) + \
                              codeError * (PDIcode / tau1code)

                    oldCodeNco = codeNco

                    oldCodeError = codeError

                    codeFreq = settings.codeFreqBasis - codeNco

                    codeFreq_[loopCnt] = codeFreq

                    # Record various measures to show in postprocessing ----------------------
                    # Record sample number (based on 8bit samples)
                    absoluteSample[loopCnt] = fid.tell()

                    dllDiscr[loopCnt] = codeError

                    dllDiscrFilt[loopCnt] = codeNco

                    pllDiscr[loopCnt] = carrError

                    pllDiscrFilt[loopCnt] = carrNco

                    I_E_[loopCnt] = I_E

                    I_P_[loopCnt] = I_P

                    I_L_[loopCnt] = I_L

                    Q_E_[loopCnt] = Q_E

                    Q_P_[loopCnt] = Q_P

                    Q_L_[loopCnt] = Q_L

                # If we got so far, this means that the tracking was successful
                # Now we only copy status, but it can be update by a lock detector
                # if implemented
                status = channel[channelNr].status
                rec.append((status, absoluteSample, codeFreq_, carrFreq_,
                            I_P_, I_E_, I_L_, Q_E_, Q_P_, Q_L_,
                            dllDiscr, dllDiscrFilt, pllDiscr, pllDiscrFilt, PRN))

        trackResults = np.rec.fromrecords(rec,
                                          dtype=[('status', 'S1'), ('absoluteSample', 'object'), ('codeFreq', 'object'),
                                                 ('carrFreq', 'object'), ('I_P', 'object'), ('I_E', 'object'),
                                                 ('I_L', 'object'),
                                                 ('Q_E', 'object'), ('Q_P', 'object'), ('Q_L', 'object'),
                                                 ('dllDiscr', 'object'),
                                                 ('dllDiscrFilt', 'object'), ('pllDiscr', 'object'),
                                                 ('pllDiscrFilt', 'object'),
                                                 ('PRN', 'int64')])
        self._results = trackResults
        return

    def plot(self):
        import matplotlib as mpl

        # %% configure matplotlib
        mpl.rcdefaults()
        # mpl.rcParams['font.sans-serif']
        # mpl.rcParams['font.family'] = 'serif'
        mpl.rc('savefig', bbox='tight', transparent=False, format='png')
        mpl.rc('axes', grid=True, linewidth=1.5, axisbelow=True)
        mpl.rc('lines', linewidth=1.5, solid_joinstyle='bevel')
        mpl.rc('figure', figsize=[8, 6], dpi=120)
        mpl.rc('text', usetex=True)
        mpl.rc('font', family='serif', serif='Computer Modern Roman', size=8)
        mpl.rc('mathtext', fontset='cm')

        # mpl.rc('font', size=16)
        # mpl.rc('text.latex', preamble=r'\usepackage{cmbright}')

        # ./plotTracking.m

        trackResults = self._results
        settings = self._settings
        channelList = range(settings.numberOfChannels)

        import matplotlib as mpl
        import matplotlib.gridspec as gs
        import matplotlib.pyplot as plt

        # %% configure matplotlib
        mpl.rcdefaults()
        # mpl.rcParams['font.sans-serif']
        # mpl.rcParams['font.family'] = 'serif'
        mpl.rc('savefig', bbox='tight', transparent=False, format='png')
        mpl.rc('axes', grid=True, linewidth=1.5, axisbelow=True)
        mpl.rc('lines', linewidth=1.5, solid_joinstyle='bevel')
        mpl.rc('figure', figsize=[8, 6], dpi=120)
        mpl.rc('text', usetex=True)
        mpl.rc('font', family='serif', serif='Computer Modern Roman', size=8)
        mpl.rc('mathtext', fontset='cm')

        # mpl.rc('font', size=16)
        # mpl.rc('text.latex', preamble=r'\usepackage{cmbright}')

        # ./plotTracking.m

        # This function plots the tracking results for the given channel list.

        # plotTracking(channelList, trackResults, settings)

        #   Inputs:
        #       channelList     - list of channels to be plotted.
        #       trackResults    - tracking results from the tracking function.
        #       settings        - receiver settings.

        # Protection - if the list contains incorrect channel numbers
        channelList = np.intersect1d(channelList, range(settings.numberOfChannels))

        # === For all listed channels ==============================================
        for channelNr in channelList:
            # Select (or create) and clear the figure ================================
            # The number 200 is added just for more convenient handling of the open
            # figure windows, when many figures are closed and reopened.
            # Figures drawn or opened by the user, will not be "overwritten" by
            # this function.
            f = plt.figure(channelNr + 200)
            f.set_label('Channel ' + str(channelNr) +
                        ' (PRN ' + str(trackResults[channelNr].PRN) + ') results')
            # Draw axes ==============================================================
            # Row 1
            spec = gs.GridSpec(3, 3)
            h11 = plt.subplot(spec[0, 0])

            h12 = plt.subplot(spec[0, 1:])

            h21 = plt.subplot(spec[1, 0])

            h22 = plt.subplot(spec[1, 1:])

            h31 = plt.subplot(spec[2, 0])

            h32 = plt.subplot(spec[2, 1])

            h33 = plt.subplot(spec[2, 2])

            # Plot all figures =======================================================
            timeAxisInSeconds = np.arange(settings.msToProcess) / 1000.0

            h11.plot(trackResults[channelNr].I_P, trackResults[channelNr].Q_P, '.')
            h11.grid()
            h11.axis('equal')
            h11.set(title='Discrete-Time Scatter Plot', xlabel='I prompt', ylabel='Q prompt')
            h12.plot(timeAxisInSeconds, trackResults[channelNr].I_P)
            h12.grid()
            h12.set(title='Bits of the navigation message', xlabel='Time (s)')
            h12.axis('tight')
            h21.plot(timeAxisInSeconds, trackResults[channelNr].pllDiscr, 'r')
            h21.grid()
            h21.axis('tight')
            h21.set(xlabel='Time (s)', ylabel='Amplitude', title='Raw PLL discriminator')
            h22.plot(timeAxisInSeconds,
                     np.sqrt(trackResults[channelNr].I_E ** 2 + trackResults[channelNr].Q_E ** 2).T,
                     timeAxisInSeconds,
                     np.sqrt(trackResults[channelNr].I_P ** 2 + trackResults[channelNr].Q_P ** 2).T,
                     timeAxisInSeconds,
                     np.sqrt(trackResults[channelNr].I_L ** 2 + trackResults[channelNr].Q_L ** 2).T, '-*')
            h22.grid()
            h22.set(title='Correlation results', xlabel='Time (s)')
            h22.axis('tight')
            h22.legend(['$\sqrt{I_{E}^2 + Q_{E}^2}$', '$\sqrt{I_{P}^2 + Q_{P}^2}$',
                        '$\sqrt{I_{L}^2 + Q_{L}^2}$'])

            h31.plot(timeAxisInSeconds, trackResults[channelNr].pllDiscrFilt, 'b')
            h31.grid()
            h31.axis('tight')
            h31.set(xlabel='Time (s)',
                    ylabel='Amplitude',
                    title='Filtered PLL discriminator')
            h32.plot(timeAxisInSeconds, trackResults[channelNr].dllDiscr, 'r')
            h32.grid()
            h32.axis('tight')
            h32.set(xlabel='Time (s)',
                    ylabel='Amplitude',
                    title='Raw DLL discriminator')
            h33.plot(timeAxisInSeconds, trackResults[channelNr].dllDiscrFilt, 'b')
            h33.grid()
            h33.axis('tight')
            h33.set(xlabel='Time (s)',
                    ylabel='Amplitude',
                    title='Filtered DLL discriminator')
            f.show()
