# Getting Started

This is a GNSS Software Defined Radio (SDR) implemented in python, based on SoftGNSS 3.0 developed by Darius Plausinaitis and Dennis M. Akos in Matlab.

# System requirements

* python 3
* matplotlib
* scipy
* numpy

# Installation

Coming soon!

# Running the GNSS SDR

1. run ```python main.py -f=your\file\path -s=(your sampling frequency in Hz) -d=(numpy dtypes, float32 .etc.)```
2. if you want to modify other settings like skipNumofSamples or set file source to be real signal not complex signal (default complex) check "initialize.py".


# Updates
* I/Q file support added
* fileType added to settings, **all datetypes** in **numpy.dtype** are supported 
* little bug solved
 
