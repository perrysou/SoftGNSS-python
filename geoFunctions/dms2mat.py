import numpy as np


# dms2mat.m
def dms2mat(dmsInput=None, n=None, *args, **kwargs):
    # DMS2MAT  Splits a real a = dd*100 + mm + s/100 into[dd mm s.ssss]
    #         where n specifies the power of 10, to which the resulting seconds
    #         of the output should be rounded. E.g.: if a result is 23.823476
    #         seconds, and n = -3, then the output will be 23.823.

    # Written by Kai Borre
    # January 7, 2007
    # Updated by Darius Plausinaitis

    neg_arg = False

    if dmsInput < 0:
        # Only positive numbers should be used while spliting into deg/min/sec
        dmsInput = -dmsInput

        neg_arg = True

    ### Split degrees minutes and seconds
    int_deg = np.floor(dmsInput / 100)

    mm = np.floor(dmsInput - 100 * int_deg)

    # we assume n<7; hence #2.10f is sufficient to hold ssdec
    ssdec = '%2.10f' % (dmsInput - 100 * int_deg - mm) * 100

    ### Check for overflow
    if ssdec == 60.0:
        mm = mm + 1

        ssdec = 0.0

    if mm == 60.0:
        int_deg = int_deg + 1

        mm = 0.0

    ### Corect the sign
    if neg_arg:
        int_deg = -int_deg

    ### Compose the output
    matOutput = []
    matOutput[0] = int_deg

    matOutput[1] = mm

    matOutput[2] = float(ssdec[0:- n + 3])

    return matOutput
################### end dms2mat.m ################
