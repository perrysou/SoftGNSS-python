import numpy as np


class AcqResults(object):
    def __init__(self):
        self.carrFreq = 0.0
        self.codePhase = 0.0
        self.peakMetric = 0.0


# generateCAcode.m
def generateCAcode(PRN=None, *args, **kwargs):
    # generateCAcode.m generates one of the 32 GPS satellite C/A codes.

    # CAcode = generateCAcode(PRN)

    #   Inputs:
    #       PRN         - PRN number of the sequence.

    #   Outputs:
    #       CAcode      - a vector containing the desired C/A code sequence
    #                   (chips).

    # --------------------------------------------------------------------------
    #                           SoftGNSS v3.0
    #
    # Copyright (C) Darius Plausinaitis
    # Written by Darius Plausinaitis
    # Based on Dennis M. Akos, Peter Rinder and Nicolaj Bertelsen
    # --------------------------------------------------------------------------
    # This program is free software; you can redistribute it and/or
    # modify it under the terms of the GNU General Public License
    # as published by the Free Software Foundation; either version 2
    # of the License, or (at your option) any later version.

    # This program is distributed in the hope that it will be useful,
    # but WITHOUT ANY WARRANTY; without even the implied warranty of
    # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    # GNU General Public License for more details.

    # You should have received a copy of the GNU General Public License
    # along with this program; if not, write to the Free Software
    # Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
    # USA.
    # --------------------------------------------------------------------------

    # CVS record:
    # $Id: generateCAcode.m,v 1.1.2.5 2006/08/14 11:38:22 dpl Exp $

    # --- Make the code shift array. The shift depends on the PRN number -------
    # The g2s vector holds the appropriate shift of the g2 code to generate
    # the C/A code (ex. for SV#19 - use a G2 shift of g2s(19) = 471)
    g2s = [5, 6, 7, 8, 17, 18, 139, 140, 141, 251,
           252, 254, 255, 256, 257, 258, 469, 470, 471, 472,
           473, 474, 509, 512, 513, 514, 515, 516, 859, 860,
           861, 862,
           145, 175, 52, 21, 237, 235, 886, 657, 634, 762, 355, 1012, 176, 603, 130, 359, 595, 68, 386]
    # generateCAcode.m:42
    # --- Pick right shift for the given PRN number ----------------------------
    g2shift = g2s[PRN]
    # generateCAcode.m:53
    # --- Generate G1 code -----------------------------------------------------

    # --- Initialize g1 output to speed up the function ---
    g1 = np.zeros(1023)
    # generateCAcode.m:58
    # --- Load shift register ---
    reg = -1 * np.ones(10)
    # generateCAcode.m:60
    # --- Generate all G1 signal chips based on the G1 feedback polynomial -----
    for i in range(1023):
        g1[i] = reg[-1]
        # generateCAcode.m:64
        saveBit = reg[2] * reg[9]
        # generateCAcode.m:65
        reg[1:] = reg[:-1]
        # generateCAcode.m:66
        reg[0] = saveBit
    # generateCAcode.m:67

    # --- Generate G2 code -----------------------------------------------------

    # --- Initialize g2 output to speed up the function ---
    g2 = np.zeros(1023)
    # generateCAcode.m:73
    # --- Load shift register ---
    reg = -1 * np.ones(10)
    # generateCAcode.m:75
    # --- Generate all G2 signal chips based on the G2 feedback polynomial -----
    for i in range(1023):
        g2[i] = reg[-1]
        # generateCAcode.m:79
        saveBit = reg[1] * reg[2] * reg[5] * reg[7] * reg[8] * reg[9]
        # generateCAcode.m:80
        reg[1:] = reg[:-1]
        # generateCAcode.m:81
        reg[0] = saveBit
    # generateCAcode.m:82

    # --- Shift G2 code --------------------------------------------------------
    # The idea: g2 = concatenate[ g2_right_part, g2_left_part ];
    g2 = np.r_[g2[1023 - g2shift:], g2[:1023 - g2shift]]
    # generateCAcode.m:87
    # --- Form single sample C/A code by multiplying G1 and G2 -----------------
    CAcode = -g1 * g2
    return CAcode


# generateCAcode.m:90


# makeCaTable.m
def makeCaTable(settings=None, *args, **kwargs):
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

    # --------------------------------------------------------------------------
    #                           SoftGNSS v3.0
    #
    # Copyright (C) Darius Plausinaitis
    # Written by Darius Plausinaitis
    # Based on Peter Rinder and Nicolaj Bertelsen
    # --------------------------------------------------------------------------
    # This program is free software; you can redistribute it and/or
    # modify it under the terms of the GNU General Public License
    # as published by the Free Software Foundation; either version 2
    # of the License, or (at your option) any later version.

    # This program is distributed in the hope that it will be useful,
    # but WITHOUT ANY WARRANTY; without even the implied warranty of
    # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    # GNU General Public License for more details.

    # You should have received a copy of the GNU General Public License
    # along with this program; if not, write to the Free Software
    # Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
    # USA.
    # --------------------------------------------------------------------------

    # CVS record:
    # $Id: makeCaTable.m,v 1.1.2.6 2006/08/14 11:38:22 dpl Exp $

    # --- Find number of samples per spreading code ----------------------------
    samplesPerCode = np.long(np.round(settings.samplingFreq / (settings.codeFreqBasis / settings.codeLength)))
    # makeCaTable.m:43
    # --- Prepare the output matrix to speed up function -----------------------
    caCodesTable = np.zeros((32, samplesPerCode))
    # makeCaTable.m:47
    # --- Find time constants --------------------------------------------------
    ts = 1.0 / settings.samplingFreq
    # makeCaTable.m:50

    tc = 1.0 / settings.codeFreqBasis
    # makeCaTable.m:51

    # === For all satellite PRN-s ...
    for PRN in range(32):
        # --- Generate CA code for given PRN -----------------------------------
        caCode = generateCAcode(PRN)
        # makeCaTable.m:56
        # --- Make index array to read C/A code values -------------------------
        # The length of the index array depends on the sampling frequency -
        # number of samples per millisecond (because one C/A code period is one
        # millisecond).
        codeValueIndex = np.ceil(ts * np.arange(1, samplesPerCode + 1) / tc) - 1
        codeValueIndex = np.longlong(codeValueIndex)
        # makeCaTable.m:64
        codeValueIndex[-1] = 1022
        # makeCaTable.m:67
        # The "upsampled" code is made by selecting values form the CA code
        # chip array (caCode) for the time instances of each sample.
        caCodesTable[PRN] = caCode[codeValueIndex]
    return caCodesTable


# makeCaTable.m:72

# ./acquisition.m
def acquisition(longSignal=None, settings=None, *args, **kwargs):
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

    # --------------------------------------------------------------------------
    #                           SoftGNSS v3.0
    #
    # Copyright (C) Darius Plausinaitis and Dennis M. Akos
    # Written by Darius Plausinaitis and Dennis M. Akos
    # Based on Peter Rinder and Nicolaj Bertelsen
    # --------------------------------------------------------------------------
    # This program is free software; you can redistribute it and/or
    # modify it under the terms of the GNU General Public License
    # as published by the Free Software Foundation; either version 2
    # of the License, or (at your option) any later version.

    # This program is distributed in the hope that it will be useful,
    # but WITHOUT ANY WARRANTY; without even the implied warranty of
    # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    # GNU General Public License for more details.

    # You should have received a copy of the GNU General Public License
    # along with this program; if not, write to the Free Software
    # Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
    # USA.
    # --------------------------------------------------------------------------

    # CVS record:
    # $Id: acquisition.m,v 1.1.2.12 2006/08/14 12:08:03 dpl Exp $

    ## Initialization =========================================================

    # Find number of samples per spreading code
    samplesPerCode = np.long(np.round(settings.samplingFreq / (settings.codeFreqBasis / settings.codeLength)))
    # ./acquisition.m:50
    # Create two 1msec vectors of data to correlate with and one with zero DC
    signal1 = longSignal[0:samplesPerCode]
    # ./acquisition.m:54
    signal2 = longSignal[samplesPerCode:2 * samplesPerCode]
    # ./acquisition.m:55
    signal0DC = longSignal - longSignal.mean()
    # ./acquisition.m:57
    # Find sampling period
    ts = 1.0 / settings.samplingFreq
    # ./acquisition.m:60
    # Find phase points of the local carrier wave
    phasePoints = np.arange(samplesPerCode) * 2 * np.pi * ts
    # ./acquisition.m:63
    # Number of the frequency bins for the given acquisition band (500Hz steps)
    numberOfFrqBins = np.int(np.round(settings.acqSearchBand * 2) + 1)
    # ./acquisition.m:66
    # Generate all C/A codes and sample them according to the sampling freq.
    caCodesTable = makeCaTable(settings)
    # ./acquisition.m:69
    # --- Initialize arrays to speed up the code -------------------------------
    # Search results of all frequency bins and code shifts (for one satellite)
    results = np.zeros((numberOfFrqBins, samplesPerCode))
    # ./acquisition.m:74
    # Carrier frequencies of the frequency bins
    frqBins = np.zeros(numberOfFrqBins)
    # ./acquisition.m:77
    # --- Initialize acqResults ------------------------------------------------
    # Carrier frequencies of detected signals
    carrFreq = np.zeros(32)
    # ./acquisition.m:82
    # C/A code phases of detected signals
    codePhase_ = np.zeros(32)
    # ./acquisition.m:84
    # Correlation peak ratios of the detected signals
    peakMetric = np.zeros(32)
    # ./acquisition.m:86
    print '('
    # Perform search for all listed PRN numbers ...
    for PRN in range(len(settings.acqSatelliteList)):
        ## Correlate signals ======================================================   
        # --- Perform DFT of C/A code ------------------------------------------
        caCodeFreqDom = np.fft.fft(caCodesTable[PRN, :]).conj()
        # ./acquisition.m:95
        for frqBinIndex in range(numberOfFrqBins):
            # --- Generate carrier wave frequency grid (0.5kHz step) -----------
            frqBins[frqBinIndex] = settings.IF - \
                                   settings.acqSearchBand / 2 * 1000 + \
                                   500.0 * frqBinIndex
            # ./acquisition.m:101
            sinCarr = np.sin(frqBins[frqBinIndex] * phasePoints)
            # ./acquisition.m:106
            cosCarr = np.cos(frqBins[frqBinIndex] * phasePoints)
            # ./acquisition.m:107
            I1 = sinCarr * signal1
            # ./acquisition.m:110
            Q1 = cosCarr * signal1
            # ./acquisition.m:111
            I2 = sinCarr * signal2
            # ./acquisition.m:112
            Q2 = cosCarr * signal2
            # ./acquisition.m:113
            IQfreqDom1 = np.fft.fft(I1 + 1j * Q1)
            # ./acquisition.m:116
            IQfreqDom2 = np.fft.fft(I2 + 1j * Q2)
            # ./acquisition.m:117
            # domain)
            convCodeIQ1 = IQfreqDom1 * caCodeFreqDom
            # ./acquisition.m:121
            convCodeIQ2 = IQfreqDom2 * caCodeFreqDom
            # ./acquisition.m:122
            acqRes1 = abs(np.fft.ifft(convCodeIQ1)) ** 2
            # ./acquisition.m:125
            acqRes2 = abs(np.fft.ifft(convCodeIQ2)) ** 2
            # ./acquisition.m:126
            # "blend" 1st and 2nd msec but will correct data bit issues
            if acqRes1.max() > acqRes2.max():
                results[frqBinIndex, :] = acqRes1
            # ./acquisition.m:131
            else:
                results[frqBinIndex, :] = acqRes2
        # ./acquisition.m:133
        ## Look for correlation peaks in the results ==============================
        # Find the highest peak and compare it to the second highest peak
        # The second peak is chosen not closer than 1 chip to the highest peak
        # --- Find the correlation peak and the carrier frequency --------------
        peakSize = results.max(1).max()
        frequencyBinIndex = results.max(1).argmax()
        # ./acquisition.m:143
        peakSize = results.max(0).max()
        codePhase = results.max(0).argmax()
        # ./acquisition.m:146
        samplesPerCodeChip = long(round(settings.samplingFreq / settings.codeFreqBasis))
        # ./acquisition.m:149
        excludeRangeIndex1 = codePhase - samplesPerCodeChip
        # ./acquisition.m:150
        excludeRangeIndex2 = codePhase + samplesPerCodeChip
        # ./acquisition.m:151
        # boundaries
        if excludeRangeIndex1 <= 0:
            codePhaseRange = np.r_[excludeRangeIndex2:samplesPerCode + excludeRangeIndex1 + 1]
        # ./acquisition.m:156
        elif excludeRangeIndex2 >= samplesPerCode - 1:
            codePhaseRange = np.r_[excludeRangeIndex2 - samplesPerCode:excludeRangeIndex1]
        # ./acquisition.m:160
        else:
            codePhaseRange = np.r_[0:excludeRangeIndex1 + 1, excludeRangeIndex2:samplesPerCode]
        # ./acquisition.m:163
        # --- Find the second highest correlation peak in the same freq. bin ---
        secondPeakSize = results[frequencyBinIndex, codePhaseRange].max()
        # ./acquisition.m:168
        peakMetric[PRN] = peakSize / secondPeakSize
        # ./acquisition.m:171
        if (peakSize / secondPeakSize) > settings.acqThreshold:
            ## Fine resolution frequency search =======================================
            # --- Indicate PRN number of the detected signal -------------------
            print '%02d ' % (PRN + 1)
            caCode = generateCAcode(PRN)
            # ./acquisition.m:182
            codeValueIndex = np.floor(ts * np.arange(1, 10 * samplesPerCode + 1) / (1.0 / settings.codeFreqBasis))
            # ./acquisition.m:184
            longCaCode = caCode[np.longlong(codeValueIndex % 1023)]
            # ./acquisition.m:187
            # (Using detected C/A code phase)
            xCarrier = signal0DC[codePhase:codePhase + 10 * samplesPerCode] * longCaCode
            # ./acquisition.m:191
            fftNumPts = 8 * 2 ** (np.ceil(np.log2(len(xCarrier))))
            # ./acquisition.m:196
            # associated carrier frequency
            fftxc = np.abs(np.fft.fft(xCarrier, np.long(fftNumPts)))
            # ./acquisition.m:200
            uniqFftPts = np.long(np.ceil((fftNumPts + 1) / 2.0))
            # ./acquisition.m:202
            fftMax = fftxc[4:uniqFftPts - 5].max()
            fftMaxIndex = fftxc[4:uniqFftPts - 5].argmax()
            # ./acquisition.m:203
            fftFreqBins = np.arange(uniqFftPts) * settings.samplingFreq / fftNumPts
            # ./acquisition.m:205
            carrFreq[PRN] = fftFreqBins[fftMaxIndex]
            # ./acquisition.m:208
            codePhase_[PRN] = codePhase
        # ./acquisition.m:209
        else:
            # --- No signal with this PRN --------------------------------------
            print '. '

    # === Acquisition is over ==================================================
    print ')\n'
    acqResults = np.core.records.fromarrays([carrFreq, codePhase_, peakMetric],
                                            names='carrFreq,codePhase,peakMetric')
    return acqResults


if __name__ == '__main__':
    acquisition()
