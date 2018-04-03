# ./initSettings.m

# Functions initializes and saves settings. Settings can be edited inside of
# the function, updated from the command line or updated using a dedicated
# GUI - "setSettings".

# All settings are described inside function code.

# settings = initSettings()

#   Inputs: none

#   Outputs:
#       settings     - Receiver settings (a structure).

class TruePosition(object):
    def __init__(self):
        self.E = None
        self.N = None
        self.U = None


class Settings(object):
    def __init__(self):
        ## Processing settings ====================================================
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

        ## Raw signal file name and other parameter ===============================
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

        ## Acquisition settings ===================================================
        # Skips acquisition in the script postProcessing.m if set to 1
        self.skipAcquisition = False

        # List of satellites to look for. Some satellites can be excluded to speed
        # up acquisition
        self.acqSatelliteList = range(1, 33)

        # Band around IF to search for satellite signal. Depends on max Doppler
        self.acqSearchBand = 14.0

        # Threshold for the signal presence decision rule
        self.acqThreshold = 2.5

        ## Tracking loops settings ================================================
        # Code tracking loop parameters
        self.dllDampingRatio = 0.7

        self.dllNoiseBandwidth = 2.0

        self.dllCorrelatorSpacing = 0.5

        # Carrier tracking loop parameters
        self.pllDampingRatio = 0.7

        self.pllNoiseBandwidth = 25.0

        ## Navigation solution settings ===========================================

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

        ## Plot settings ==========================================================
        # Enable/disable plotting of the tracking results for each channel
        self.plotTracking = True

        # 1 - On

        ## Constants ==============================================================

        self.c = 299792458.0

        self.startOffset = 68.802
