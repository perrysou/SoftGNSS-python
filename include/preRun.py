import numpy as np
from scipy.io import loadmat


class Channel(object):
    def __init__(self):
        self.PRN = 0
        # preRun.m:46
        self.acquiredFreq = 0.0
        # preRun.m:47
        self.codePhase = 0.0
        # preRun.m:48
        self.status = '-'
        # preRun.m:50
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
    # $Id: preRun.m,v 1.8.2.20 2006/08/14 11:38:22 dpl Exp $

    ## Initialize all channels ================================================
    PRN = np.zeros(settings.numberOfChannels, dtype='int64')
    acquiredFreq = np.zeros(settings.numberOfChannels)
    codePhase = np.zeros(settings.numberOfChannels)
    status = ['-' for _ in range(settings.numberOfChannels)]
    # preRun.m:44
    # --- Copy initial data to all channels ------------------------------------
    # preRun.m:55
    ## Copy acquisition results ===============================================

    # --- Sort peaks to find strongest signals, keep the peak index information
    PRNindexes = sorted(enumerate(acqResults.peakMetric),
                        key=lambda x: x[-1], reverse=True)
    # preRun.m:60
    # --- Load information about each satellite --------------------------------
    # Maximum number of initialized channels is number of detected signals, but
    # not more as the number of channels specified in the settings.
    for ii in range(min(settings.numberOfChannels, sum(acqResults.carrFreq > 0))):
        PRN[ii] = PRNindexes[ii][0] + 1
        # preRun.m:66
        acquiredFreq[ii] = acqResults.carrFreq[PRNindexes[ii][0]]
        # preRun.m:67
        codePhase[ii] = acqResults.codePhase[PRNindexes[ii][0]]
        # preRun.m:68
        status[ii] = 'T'
    # preRun.m:71
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

    # --------------------------------------------------------------------------
    #                           SoftGNSS v3.0
    #
    # Copyright (C) Peter Rinder and Nicolaj Bertelsen
    # Written by Peter Rinder Nicolaj Bertelsen and Darius Plausinaitis
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
    # $Id: showChannelStatus.m,v 1.4.2.8 2006/08/14 11:38:22 dpl Exp $

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
