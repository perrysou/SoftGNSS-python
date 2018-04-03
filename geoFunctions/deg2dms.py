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

    if deg < 0:
        # Only positive numbers should be used while spliting into deg/min/sec
        deg = -deg

        neg_arg = True

    ### Split degrees minutes and seconds
    int_deg = np.floor(deg)

    decimal = deg - int_deg

    min_part = decimal * 60

    min_ = np.floor(min_part)

    sec_part = min_part - np.floor(min_part)

    sec = sec_part * 60

    ### Check for overflow
    if sec == 60.0:
        min_ = min_ + 1

        sec = 0.0

    if min_ == 60.0:
        int_deg = int_deg + 1

        min_ = 0.0

    ### Construct the output
    dmsOutput = int_deg * 100 + min_ + sec / 100

    ### Correct the sign
    if neg_arg:
        dmsOutput = -dmsOutput
    return dmsOutput

################### end deg2dms.m ################
