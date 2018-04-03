import numpy as np
from scipy.io import loadmat


class Channel(object):
    def __init__(self):
        self.PRN = 0

        self.acquiredFreq = 0.0

        self.codePhase = 0.0

        self.status = '-'

        # "-" - "off" - no signal to track
        # "T" - Tracking state


# preRun.m


def preRun(acqResults=None, settings=None, *args, **kwargs):
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

    ## Initialize all channels ================================================
    PRN = np.zeros(settings.numberOfChannels, dtype='int64')
    acquiredFreq = np.zeros(settings.numberOfChannels)
    codePhase = np.zeros(settings.numberOfChannels)
    status = ['-' for _ in range(settings.numberOfChannels)]

    # --- Copy initial data to all channels ------------------------------------

    ## Copy acquisition results ===============================================

    # --- Sort peaks to find strongest signals, keep the peak index information
    PRNindexes = sorted(enumerate(acqResults.peakMetric),
                        key=lambda x: x[-1], reverse=True)

    # --- Load information about each satellite --------------------------------
    # Maximum number of initialized channels is number of detected signals, but
    # not more as the number of channels specified in the settings.
    for ii in range(min(settings.numberOfChannels, sum(acqResults.carrFreq > 0))):
        PRN[ii] = PRNindexes[ii][0] + 1

        acquiredFreq[ii] = acqResults.carrFreq[PRNindexes[ii][0]]

        codePhase[ii] = acqResults.codePhase[PRNindexes[ii][0]]

        status[ii] = 'T'

    channel = np.core.records.fromarrays([PRN, acquiredFreq, codePhase, status],
                                         names='PRN,acquiredFreq,codePhase,status')
    return channel


# showChannelStatus.m


def showChannelStatus(channel=None, settings=None, *args, **kwargs):
    # Prints the status of all channels in a table.

    # showChannelStatus(channel, settings)

    #   Inputs:
    #       channel     - data for each channel. It is used to initialize and
    #                   at the processing of the signal (tracking part).
    #       settings    - receiver settings

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
    a = loadmat('../trackingResults.mat', struct_as_record=False, squeeze_me=True)
    acqResults = a['acqResults']
    settings = a['settings']
    channel = preRun(acqResults, settings)
    showChannelStatus(channel, settings)
