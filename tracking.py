import numpy as np

import acquisition


class TrackResults(object):
    def __init__(self, msToProcess):
        ## Initialize result structure ============================================
        # Channel status
        self.status = '-'

        # The absolute sample in the record of the C/A code start:
        self.absoluteSample = np.zeros(msToProcess)

        # Freq of the C/A code:
        self.codeFreq = np.Inf * np.ones(msToProcess)

        # Frequency of the tracked carrier wave:
        self.carrFreq = np.Inf * np.ones(msToProcess)

        # Outputs from the correlators (In-phase):
        self.I_P = np.zeros(msToProcess)

        self.I_E = np.zeros(msToProcess)

        self.I_L = np.zeros(msToProcess)

        # Outputs from the correlators (Quadrature-phase):
        self.Q_E = np.zeros(msToProcess)

        self.Q_P = np.zeros(msToProcess)

        self.Q_L = np.zeros(msToProcess)

        # Loop discriminators
        self.dllDiscr = np.Inf * np.ones(msToProcess)

        self.dllDiscrFilt = np.Inf * np.ones(msToProcess)

        self.pllDiscr = np.Inf * np.ones(msToProcess)

        self.pllDiscrFilt = np.Inf * np.ones(msToProcess)

        self.PRN = 0


# calcLoopCoef.m


def calcLoopCoef(LBW, zeta, k, *args, **kwargs):
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


# ./tracking.m


def tracking(fid=None, channel=None, settings=None, *args, **kwargs):
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

    ## Initialize tracking variables ==========================================

    codePeriods = settings.msToProcess

    # --- DLL variables --------------------------------------------------------
    # Define early-late offset (in chips)
    earlyLateSpc = settings.dllCorrelatorSpacing

    # Summation interval
    PDIcode = 0.001

    # Calculate filter coefficient values
    tau1code, tau2code = calcLoopCoef(settings.dllNoiseBandwidth, settings.dllDampingRatio, 1.0)

    # --- PLL variables --------------------------------------------------------
    # Summation interval
    PDIcarr = 0.001

    # Calculate filter coefficient values
    tau1carr, tau2carr = calcLoopCoef(settings.pllNoiseBandwidth, settings.pllDampingRatio, 0.25)

    #     hwb=waitbar(0,'Tracking...')

    # Initialize a temporary list of records
    rec = []
    ## Start processing channels ==============================================
    for channelNr in range(settings.numberOfChannels):
        settings.msToProcess = np.long(settings.msToProcess)
        # Initialize fields for record(structured) array of tracked results
        status = '-'

        # The absolute sample in the record of the C/A code start:
        absoluteSample = np.zeros(settings.msToProcess)

        # Freq of the C/A code:
        codeFreq_ = np.Inf * np.ones(settings.msToProcess)

        # Frequency of the tracked carrier wave:
        carrFreq_ = np.Inf * np.ones(settings.msToProcess)

        # Outputs from the correlators (In-phase):
        I_P_ = np.zeros(settings.msToProcess)

        I_E_ = np.zeros(settings.msToProcess)

        I_L_ = np.zeros(settings.msToProcess)

        # Outputs from the correlators (Quadrature-phase):
        Q_E_ = np.zeros(settings.msToProcess)

        Q_P_ = np.zeros(settings.msToProcess)

        Q_L_ = np.zeros(settings.msToProcess)

        # Loop discriminators
        dllDiscr = np.Inf * np.ones(settings.msToProcess)

        dllDiscrFilt = np.Inf * np.ones(settings.msToProcess)

        pllDiscr = np.Inf * np.ones(settings.msToProcess)

        pllDiscrFilt = np.Inf * np.ones(settings.msToProcess)

        PRN = 0

        settings.msToProcess = np.longfloat(settings.msToProcess)
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
            caCode = acquisition.generateCAcode(channel[channelNr].PRN - 1)

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
                ## GUI update -------------------------------------------------------------
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
                ## Read next block of data ------------------------------------------------            
                # Find the size of a "block" or code period in whole samples
                # Update the phasestep based on code freq (variable) and
                # sampling frequency (fixed)
                codePhaseStep = codeFreq / settings.samplingFreq

                blksize = np.ceil((settings.codeLength - remCodePhase) / codePhaseStep)
                blksize = np.long(blksize)

                # interation
                rawSignal = np.fromfile(fid, settings.dataType, blksize)
                samplesRead = len(rawSignal)

                # If did not read in enough samples, then could be out of
                # data - better exit
                if samplesRead != blksize:
                    print 'Not able to read the specified number of samples for tracking, exiting!'
                    fid.close()
                    trackResults = None
                    return trackResults, channel
                ## Set up all the code phase tracking information -------------------------
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

                ## Generate the carrier frequency to mix the signal to baseband -----------
                time = np.arange(0, blksize + 1) / settings.samplingFreq

                trigarg = carrFreq * 2.0 * np.pi * time + remCarrPhase

                remCarrPhase = trigarg[blksize] % (2 * np.pi)

                carrCos = np.cos(trigarg[0:blksize])

                carrSin = np.sin(trigarg[0:blksize])

                ## Generate the six standard accumulated values ---------------------------
                # First mix to baseband
                qBasebandSignal = carrCos * rawSignal

                iBasebandSignal = carrSin * rawSignal

                I_E = (earlyCode * iBasebandSignal).sum()

                Q_E = (earlyCode * qBasebandSignal).sum()

                I_P = (promptCode * iBasebandSignal).sum()

                Q_P = (promptCode * qBasebandSignal).sum()

                I_L = (lateCode * iBasebandSignal).sum()

                Q_L = (lateCode * qBasebandSignal).sum()

                ## Find PLL error and update carrier NCO ----------------------------------
                # Implement carrier loop discriminator (phase detector)
                carrError = np.arctan(Q_P / I_P) / 2.0 / np.pi

                carrNco = oldCarrNco + \
                          tau2carr / tau1carr * (carrError - oldCarrError) + \
                          carrError * (PDIcarr / tau1carr)

                oldCarrNco = carrNco

                oldCarrError = carrError

                carrFreq = carrFreqBasis + carrNco

                carrFreq_[loopCnt] = carrFreq

                ## Find DLL error and update code NCO -------------------------------------
                codeError = (np.sqrt(I_E * I_E + Q_E * Q_E) - np.sqrt(I_L * I_L + Q_L * Q_L)) / (
                        np.sqrt(I_E * I_E + Q_E * Q_E) + np.sqrt(I_L * I_L + Q_L * Q_L))

                codeNco = oldCodeNco + \
                          tau2code / tau1code * (codeError - oldCodeError) + \
                          codeError * (PDIcode / tau1code)

                oldCodeNco = codeNco

                oldCodeError = codeError

                codeFreq = settings.codeFreqBasis - codeNco

                codeFreq_[loopCnt] = codeFreq

                ## Record various measures to show in postprocessing ----------------------
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

    # Close the waitbar
    # close_(hwb)
    trackResults = np.rec.fromrecords(rec,
                                      dtype=[('status', 'S1'), ('absoluteSample', 'object'), ('codeFreq', 'object'),
                                             ('carrFreq', 'object'), ('I_P', 'object'), ('I_E', 'object'),
                                             ('I_L', 'object'),
                                             ('Q_E', 'object'), ('Q_P', 'object'), ('Q_L', 'object'),
                                             ('dllDiscr', 'object'),
                                             ('dllDiscrFilt', 'object'), ('pllDiscr', 'object'),
                                             ('pllDiscrFilt', 'object'),
                                             ('PRN', 'int64')])
    return trackResults, channel
