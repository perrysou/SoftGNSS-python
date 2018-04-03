import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import hamming, welch


# ./probeData.m


# @function
def probeData(*args, **kwargs):
    nargin = len(args)

    # Function plots raw data information: time domain plot, a frequency domain
    # plot and a histogram.

    # The function can be called in two ways:
    #   probeData(settings)
    # or
    #   probeData(fileName, settings)

    #   Inputs:
    #       fileName        - name of the data file. File name is read from
    #                       settings if parameter fileName is not provided.

    #       settings        - receiver settings. Type of data file, sampling
    #                       frequency and the default filename are specified
    #                       here.

    ## Check the number of arguments ==========================================
    if nargin == 1:
        settings = args[0]

        fileNameStr = settings.fileName

    elif nargin == 2:
        fileNameStr, settings = args

        if not isinstance(fileNameStr, str):
            raise TypeError('File name must be a string')
    else:
        raise Exception('Incorrect number of arguments')

    ## Generate plot of raw data ==============================================

    try:
        with open(fileNameStr, 'rb') as fid:
            # Move the starting point of processing. Can be used to start the
            # signal processing at any point in the data record (e.g. for long
            # records).
            fid.seek(settings.skipNumberOfBytes, 0)
            samplesPerCode = long(round(settings.samplingFreq / (settings.codeFreqBasis / settings.codeLength)))

            try:
                data = np.fromfile(fid,
                                   settings.dataType,
                                   10 * samplesPerCode)

            except IOError:
                # The file is too short
                print 'Could not read enough data from the data file.'
            # --- Initialization ---------------------------------------------------
            plt.figure(100)
            plt.clf()
            timeScale = np.arange(0, 0.005, 1 / settings.samplingFreq)

            plt.subplot(2, 2, 1)
            plt.plot(1000 * timeScale[1:round(samplesPerCode / 50)],
                     data[1:round(samplesPerCode / 50)])
            plt.axis('tight')
            plt.grid()
            plt.title('Time domain plot')
            plt.xlabel('Time (ms)')
            plt.ylabel('Amplitude')
            plt.subplot(2, 2, 2)
            f, Pxx = welch(data - np.mean(data),
                           settings.samplingFreq / 1000000.0,
                           hamming(16384, False),
                           16384,
                           1024,
                           16384)
            plt.semilogy(f, Pxx)
            plt.axis('tight')
            plt.grid()
            plt.title('Frequency domain plot')
            plt.xlabel('Frequency (MHz)')
            plt.ylabel('Magnitude')
            plt.show()
            plt.subplot(2, 2, 3.5)
            plt.hist(data, np.arange(- 128, 128))
            dmax = np.max(np.abs(data)) + 1

            plt.axis('tight')
            adata = plt.axis()

            plt.axis([-dmax, dmax, adata[2], adata[3]])
            plt.grid('on')
            plt.title('Histogram')
            plt.xlabel('Bin')
            plt.ylabel('Number in bin')
        # === Error while opening the data file ================================
    except IOError as e:
        print 'Unable to read file "%s": %s' % (fileNameStr, e)
