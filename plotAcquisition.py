import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.io.matlab import loadmat

# %% configure matplotlib
mpl.rcdefaults()
# mpl.rcParams['font.sans-serif']
# mpl.rcParams['font.family'] = 'serif'
mpl.rc('savefig', bbox='tight', transparent=False, format='png')
mpl.rc('axes', grid=True, linewidth=1.5, axisbelow=True)
mpl.rc('lines', linewidth=1.5, solid_joinstyle='bevel')
mpl.rc('figure', figsize=[8, 6], autolayout=False, dpi=120)
mpl.rc('text', usetex=True)
mpl.rc('font', family='serif', serif='computer modern roman', size=16)
mpl.rc('mathtext', fontset='cm')


# mpl.rc('font', size=16)
# mpl.rc('text.latex', preamble=r'\usepackage{cmbright}')

# ./plotAcquisition.m


def plotAcquisition(acqResults=None):
    # Functions plots bar plot of acquisition results (acquisition metrics). No
    # bars are shown for the satellites not included in the acquisition list (in
    # structure SETTINGS).

    # plotAcquisition(acqResults)

    #   Inputs:
    #       acqResults    - Acquisition results from function acquisition.

    ## Plot all results =======================================================
    f, hAxes = plt.subplots()
    # ./plotAcquisition.m:39
    plt.bar(range(1, 33), acqResults.peakMetric)
    plt.title('Acquisition results')
    plt.xlabel('PRN number (no bar - SV is not in the acquisition list)')
    plt.ylabel('Acquisition Metric ($1^{st}$ to $2^{nd}$ Correlation Peaks Ratio')
    oldAxis = plt.axis()
    # ./plotAcquisition.m:47
    plt.axis([0, 33, 0, oldAxis[-1]])
    plt.xticks(range(1, 33), size=12)
    # plt.minorticks_on()
    hAxes.xaxis.grid()
    ## Mark acquired signals ==================================================

    acquiredSignals = acqResults.peakMetric * (acqResults.carrFreq > 0)
    # ./plotAcquisition.m:54
    plt.bar(range(1, 33), acquiredSignals, FaceColor=(0, 0.8, 0))
    plt.legend(['Not acquired signals', 'Acquired signals'])
    plt.show()


if __name__ == '__main__':
    a = loadmat('./trackingResults.mat', struct_as_record=False, squeeze_me=True)
    acqResults = a['acqResults']
    plotAcquisition(acqResults)
