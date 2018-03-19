import numpy as np


# deg2dms.m
def deg2dms(deg=None, *args, **kwargs):
    # DEG2DMS  Conversion of degrees to degrees, minutes, and seconds.
    # The output format (dms format) is: (degrees*100 + minutes + seconds/100)

    # Written by Kai Borre
    # February 7, 2001
    # Updated by Darius Plausinaitis

    ### Save the sign for later processing
    neg_arg = False
    # deg2dms.m:10
    if deg < 0:
        # Only positive numbers should be used while spliting into deg/min/sec
        deg = -deg
        # deg2dms.m:13
        neg_arg = True
    # deg2dms.m:14

    ### Split degrees minutes and seconds
    int_deg = np.floor(deg)
    # deg2dms.m:18
    decimal = deg - int_deg
    # deg2dms.m:19
    min_part = decimal * 60
    # deg2dms.m:20
    min_ = np.floor(min_part)
    # deg2dms.m:21
    sec_part = min_part - np.floor(min_part)
    # deg2dms.m:22
    sec = sec_part * 60
    # deg2dms.m:23
    ### Check for overflow
    if sec == 60.0:
        min_ = min_ + 1
        # deg2dms.m:27
        sec = 0.0
    # deg2dms.m:28

    if min_ == 60.0:
        int_deg = int_deg + 1
        # deg2dms.m:31
        min_ = 0.0
    # deg2dms.m:32

    ### Construct the output
    dmsOutput = int_deg * 100 + min_ + sec / 100
    # deg2dms.m:36
    ### Correct the sign
    if neg_arg:
        dmsOutput = -dmsOutput
    return dmsOutput
# deg2dms.m:40

################### end deg2dms.m ################
