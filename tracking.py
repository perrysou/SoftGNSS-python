import numpy as np

import acquisition


class TrackResults(object):
    def __init__(self, msToProcess):
        ## Initialize result structure ============================================
        # Channel status
        self.status = '-'
        # ./tracking.m:48
        # The absolute sample in the record of the C/A code start:
        self.absoluteSample = np.zeros(msToProcess)
        # ./tracking.m:51
        # Freq of the C/A code:
        self.codeFreq = np.Inf * np.ones(msToProcess)
        # ./tracking.m:54
        # Frequency of the tracked carrier wave:
        self.carrFreq = np.Inf * np.ones(msToProcess)
        # ./tracking.m:57
        # Outputs from the correlators (In-phase):
        self.I_P = np.zeros(msToProcess)
        # ./tracking.m:60
        self.I_E = np.zeros(msToProcess)
        # ./tracking.m:61
        self.I_L = np.zeros(msToProcess)
        # ./tracking.m:62
        # Outputs from the correlators (Quadrature-phase):
        self.Q_E = np.zeros(msToProcess)
        # ./tracking.m:65
        self.Q_P = np.zeros(msToProcess)
        # ./tracking.m:66
        self.Q_L = np.zeros(msToProcess)
        # ./tracking.m:67
        # Loop discriminators
        self.dllDiscr = np.Inf * np.ones(msToProcess)
        # ./tracking.m:70
        self.dllDiscrFilt = np.Inf * np.ones(msToProcess)
        # ./tracking.m:71
        self.pllDiscr = np.Inf * np.ones(msToProcess)
        # ./tracking.m:72
        self.pllDiscrFilt = np.Inf * np.ones(msToProcess)
        # ./tracking.m:73
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
    # calcLoopCoef.m:41
    # solve for t1 & t2
    tau1 = k / (Wn * Wn)
    # calcLoopCoef.m:44
    tau2 = 2.0 * zeta / Wn
    # calcLoopCoef.m:45
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
    # ./tracking.m:80

    # --- DLL variables --------------------------------------------------------
    # Define early-late offset (in chips)
    earlyLateSpc = settings.dllCorrelatorSpacing
    # ./tracking.m:84
    # Summation interval
    PDIcode = 0.001
    # ./tracking.m:87
    # Calculate filter coefficient values
    tau1code, tau2code = calcLoopCoef(settings.dllNoiseBandwidth, settings.dllDampingRatio, 1.0)
    # ./tracking.m:90
    # --- PLL variables --------------------------------------------------------
    # Summation interval
    PDIcarr = 0.001
    # ./tracking.m:96
    # Calculate filter coefficient values
    tau1carr, tau2carr = calcLoopCoef(settings.pllNoiseBandwidth, settings.pllDampingRatio, 0.25)
    # ./tracking.m:99
    #     hwb=waitbar(0,'Tracking...')
    # ./tracking.m:102
    # Initialize a temporary list of records
    rec = []
    ## Start processing channels ==============================================
    for channelNr in range(settings.numberOfChannels):
        settings.msToProcess = np.long(settings.msToProcess)
        # Initialize fields for record(structured) array of tracked results
        status = '-'
        # ./tracking.m:48
        # The absolute sample in the record of the C/A code start:
        absoluteSample = np.zeros(settings.msToProcess)
        # ./tracking.m:51
        # Freq of the C/A code:
        codeFreq_ = np.Inf * np.ones(settings.msToProcess)
        # ./tracking.m:54
        # Frequency of the tracked carrier wave:
        carrFreq_ = np.Inf * np.ones(settings.msToProcess)
        # ./tracking.m:57
        # Outputs from the correlators (In-phase):
        I_P_ = np.zeros(settings.msToProcess)
        # ./tracking.m:60
        I_E_ = np.zeros(settings.msToProcess)
        # ./tracking.m:61
        I_L_ = np.zeros(settings.msToProcess)
        # ./tracking.m:62
        # Outputs from the correlators (Quadrature-phase):
        Q_E_ = np.zeros(settings.msToProcess)
        # ./tracking.m:65
        Q_P_ = np.zeros(settings.msToProcess)
        # ./tracking.m:66
        Q_L_ = np.zeros(settings.msToProcess)
        # ./tracking.m:67
        # Loop discriminators
        dllDiscr = np.Inf * np.ones(settings.msToProcess)
        # ./tracking.m:70
        dllDiscrFilt = np.Inf * np.ones(settings.msToProcess)
        # ./tracking.m:71
        pllDiscr = np.Inf * np.ones(settings.msToProcess)
        # ./tracking.m:72
        pllDiscrFilt = np.Inf * np.ones(settings.msToProcess)
        # ./tracking.m:73
        PRN = 0
        # ./tracking.m:76
        settings.msToProcess = np.longfloat(settings.msToProcess)
        # Only process if PRN is non zero (acquisition was successful)
        if channel[channelNr].PRN != 0:
            # Save additional information - each channel's tracked PRN
            PRN = channel[channelNr].PRN
            # ./tracking.m:110
            # signal processing at any point in the data record (e.g. for long
            # records). In addition skip through that data file to start at the
            # appropriate sample (corresponding to code phase). Assumes sample
            # type is schar (or 1 byte per sample)
            fid.seek(settings.skipNumberOfBytes + channel[channelNr].codePhase, 0)
            # Here PRN is the actual satellite ID instead of the 0-based index
            caCode = acquisition.generateCAcode(channel[channelNr].PRN - 1)
            # ./tracking.m:123
            caCode = np.r_[caCode[-1], caCode, caCode[0]]
            # ./tracking.m:125
            # define initial code frequency basis of NCO
            codeFreq = settings.codeFreqBasis
            # ./tracking.m:130
            remCodePhase = 0.0
            # ./tracking.m:132
            carrFreq = channel[channelNr].acquiredFreq
            # ./tracking.m:134
            carrFreqBasis = channel[channelNr].acquiredFreq
            # ./tracking.m:135
            remCarrPhase = 0.0
            # ./tracking.m:137
            oldCodeNco = 0.0
            # ./tracking.m:140
            oldCodeError = 0.0
            # ./tracking.m:141
            oldCarrNco = 0.0
            # ./tracking.m:144
            oldCarrError = 0.0
            # ./tracking.m:145
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
                # ./tracking.m:176
                blksize = np.ceil((settings.codeLength - remCodePhase) / codePhaseStep)
                blksize = np.long(blksize)
                # ./tracking.m:178
                # interation
                rawSignal = np.fromfile(fid, settings.dataType, blksize)
                samplesRead = len(rawSignal)
                # ./tracking.m:182
                # ./tracking.m:184
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
                # ./tracking.m:196
                tcode2 = np.ceil(tcode)
                # ./tracking.m:199
                earlyCode = caCode[np.longlong(tcode2)]
                # ./tracking.m:200
                tcode = np.linspace(remCodePhase + earlyLateSpc,
                                    blksize * codePhaseStep + remCodePhase + earlyLateSpc,
                                    blksize, endpoint=False)
                # ./tracking.m:203
                tcode2 = np.ceil(tcode)
                # ./tracking.m:206
                lateCode = caCode[np.longlong(tcode2)]
                # ./tracking.m:207
                tcode = np.linspace(remCodePhase,
                                    blksize * codePhaseStep + remCodePhase,
                                    blksize, endpoint=False)
                # ./tracking.m:210
                tcode2 = np.ceil(tcode)
                # ./tracking.m:213
                promptCode = caCode[np.longlong(tcode2)]
                # ./tracking.m:214
                remCodePhase = tcode[blksize - 1] + codePhaseStep - 1023.0
                # ./tracking.m:216
                ## Generate the carrier frequency to mix the signal to baseband -----------
                time = np.arange(0, blksize + 1) / settings.samplingFreq
                # ./tracking.m:219
                trigarg = carrFreq * 2.0 * np.pi * time + remCarrPhase
                # ./tracking.m:222
                remCarrPhase = trigarg[blksize] % (2 * np.pi)
                # ./tracking.m:223
                carrCos = np.cos(trigarg[0:blksize])
                # ./tracking.m:226
                carrSin = np.sin(trigarg[0:blksize])
                # ./tracking.m:227
                ## Generate the six standard accumulated values ---------------------------
                # First mix to baseband
                qBasebandSignal = carrCos * rawSignal
                # ./tracking.m:231
                iBasebandSignal = carrSin * rawSignal
                # ./tracking.m:232
                I_E = (earlyCode * iBasebandSignal).sum()
                # ./tracking.m:235
                Q_E = (earlyCode * qBasebandSignal).sum()
                # ./tracking.m:236
                I_P = (promptCode * iBasebandSignal).sum()
                # ./tracking.m:237
                Q_P = (promptCode * qBasebandSignal).sum()
                # ./tracking.m:238
                I_L = (lateCode * iBasebandSignal).sum()
                # ./tracking.m:239
                Q_L = (lateCode * qBasebandSignal).sum()
                # ./tracking.m:240
                ## Find PLL error and update carrier NCO ----------------------------------
                # Implement carrier loop discriminator (phase detector)
                carrError = np.arctan(Q_P / I_P) / 2.0 / np.pi
                # ./tracking.m:245
                carrNco = oldCarrNco + \
                          tau2carr / tau1carr * (carrError - oldCarrError) + \
                          carrError * (PDIcarr / tau1carr)
                # ./tracking.m:248
                oldCarrNco = carrNco
                # ./tracking.m:250
                oldCarrError = carrError
                # ./tracking.m:251
                carrFreq = carrFreqBasis + carrNco
                # ./tracking.m:254
                carrFreq_[loopCnt] = carrFreq
                # ./tracking.m:256
                ## Find DLL error and update code NCO -------------------------------------
                codeError = (np.sqrt(I_E * I_E + Q_E * Q_E) - np.sqrt(I_L * I_L + Q_L * Q_L)) / (
                        np.sqrt(I_E * I_E + Q_E * Q_E) + np.sqrt(I_L * I_L + Q_L * Q_L))
                # ./tracking.m:259
                codeNco = oldCodeNco + \
                          tau2code / tau1code * (codeError - oldCodeError) + \
                          codeError * (PDIcode / tau1code)
                # ./tracking.m:263
                oldCodeNco = codeNco
                # ./tracking.m:265
                oldCodeError = codeError
                # ./tracking.m:266
                codeFreq = settings.codeFreqBasis - codeNco
                # ./tracking.m:269
                codeFreq_[loopCnt] = codeFreq
                # ./tracking.m:271
                ## Record various measures to show in postprocessing ----------------------
                # Record sample number (based on 8bit samples)
                absoluteSample[loopCnt] = fid.tell()
                # ./tracking.m:275
                dllDiscr[loopCnt] = codeError
                # ./tracking.m:277
                dllDiscrFilt[loopCnt] = codeNco
                # ./tracking.m:278
                pllDiscr[loopCnt] = carrError
                # ./tracking.m:279
                pllDiscrFilt[loopCnt] = carrNco
                # ./tracking.m:280
                I_E_[loopCnt] = I_E
                # ./tracking.m:282
                I_P_[loopCnt] = I_P
                # ./tracking.m:283
                I_L_[loopCnt] = I_L
                # ./tracking.m:284
                Q_E_[loopCnt] = Q_E
                # ./tracking.m:285
                Q_P_[loopCnt] = Q_P
                # ./tracking.m:286
                Q_L_[loopCnt] = Q_L
            # ./tracking.m:287
            # If we got so far, this means that the tracking was successful
            # Now we only copy status, but it can be update by a lock detector
            # if implemented
            status = channel[channelNr].status
            rec.append((status, absoluteSample, codeFreq_, carrFreq_,
                        I_P_, I_E_, I_L_, Q_E_, Q_P_, Q_L_,
                        dllDiscr, dllDiscrFilt, pllDiscr, pllDiscrFilt, PRN))
    # ./tracking.m:293

    # Close the waitbar
    # close_(hwb)
    trackResults = np.rec.fromrecords(rec,
                               dtype=[('status', 'S1'), ('absoluteSample', 'object'), ('codeFreq', 'object'),
                                       ('carrFreq', 'object'), ('I_P', 'object'), ('I_E', 'object'), ('I_L', 'object'),
                                       ('Q_E', 'object'), ('Q_P', 'object'), ('Q_L', 'object'), ('dllDiscr', 'object'),
                                       ('dllDiscrFilt', 'object'), ('pllDiscr', 'object'), ('pllDiscrFilt', 'object'),
                                       ('PRN', 'int64')])
    return trackResults, channel
