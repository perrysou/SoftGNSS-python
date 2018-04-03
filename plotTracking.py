import matplotlib as mpl
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
import numpy as np
from scipy.io.matlab import loadmat

# %% configure matplotlib
mpl.rcdefaults()
# mpl.rcParams['font.sans-serif']
# mpl.rcParams['font.family'] = 'serif'
mpl.rc('savefig', bbox='tight', transparent=False, format='png')
mpl.rc('axes', grid=True, linewidth=1.5, axisbelow=True)
mpl.rc('lines', linewidth=1.5, solid_joinstyle='bevel')
mpl.rc('figure', figsize=[8, 6], dpi=120)
mpl.rc('text', usetex=True)
mpl.rc('font', family='serif', serif='computer modern roman', size=8)
mpl.rc('mathtext', fontset='cm')


# mpl.rc('font', size=16)
# mpl.rc('text.latex', preamble=r'\usepackage{cmbright}')

# ./plotTracking.m


def plotTracking(channelList=None, trackResults=None, settings=None, *args, **kwargs):
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
        ## Select (or create) and clear the figure ================================
        # The number 200 is added just for more convenient handling of the open
        # figure windows, when many figures are closed and reopened.
        # Figures drawn or opened by the user, will not be "overwritten" by
        # this function.
        f = plt.figure(channelNr + 200)
        f.set_label('Channel ' + str(channelNr) +
                    ' (PRN ' + str(trackResults[channelNr].PRN) + ') results')
        ## Draw axes ==============================================================
        # Row 1
        spec = gs.GridSpec(3, 3)
        h11 = plt.subplot(spec[0, 0])

        h12 = plt.subplot(spec[0, 1:])

        h21 = plt.subplot(spec[1, 0])

        h22 = plt.subplot(spec[1, 1:])

        h31 = plt.subplot(spec[2, 0])

        h32 = plt.subplot(spec[2, 1])

        h33 = plt.subplot(spec[2, 2])

        ## Plot all figures =======================================================
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


if __name__ == '__main__':
    a = loadmat('./trackingResults.mat', struct_as_record=False, squeeze_me=True)
    trackResults = a['trackResults']
    settings = a['settings']
    channelList = range(a['channel'].size)
    plotTracking(channelList, trackResults, settings)
