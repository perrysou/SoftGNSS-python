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
    # dms2mat.m:11
    if dmsInput < 0:
        # Only positive numbers should be used while spliting into deg/min/sec
        dmsInput = -dmsInput
        # dms2mat.m:14
        neg_arg = True
    # dms2mat.m:15

    ### Split degrees minutes and seconds
    int_deg = np.floor(dmsInput / 100)
    # dms2mat.m:19
    mm = np.floor(dmsInput - 100 * int_deg)
    # dms2mat.m:20
    # we assume n<7; hence #2.10f is sufficient to hold ssdec
    ssdec = '%2.10f' % (dmsInput - 100 * int_deg - mm) * 100
    # dms2mat.m:22
    ### Check for overflow
    if ssdec == 60.0:
        mm = mm + 1
        # dms2mat.m:26
        ssdec = 0.0
    # dms2mat.m:27

    if mm == 60.0:
        int_deg = int_deg + 1
        # dms2mat.m:30
        mm = 0.0
    # dms2mat.m:31

    ### Corect the sign
    if neg_arg:
        int_deg = -int_deg
    # dms2mat.m:36

    ### Compose the output
    matOutput = []
    matOutput[0] = int_deg
    # dms2mat.m:40
    matOutput[1] = mm
    # dms2mat.m:41
    matOutput[2] = float(ssdec[0:- n + 3])

    # dms2mat.m:42
    return matOutput
################### end dms2mat.m ################
