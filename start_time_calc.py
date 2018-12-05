#!/usr/bin/env python
# calculates the start time in seconds and MJD of a given RAW file.
# call : python start_time_calc.py filename.raw
# (code follows https://github.com/UCBerkeleySETI/breakthrough/blob/master/doc/Time-in-RAW-Files.md)

import sys
import math
import os


def read_header(fname):
    with open(fname, 'rb') as f:
        currline = str(f.read(80))

        while currline[0:3] != 'END':               # until reaching end of header
            currline = str(f.read(80))              # read new header line
            if currline.startswith('STT_IMJD='):    # read start MJD integer part
                STT_IMJD = float(currline[9:])
            if currline.startswith('STT_SMJD='):    # read start MJD seconds part
                STT_SMJD = float(currline[9:])
            if currline.startswith('PKTIDX  ='):    # read block size
                PKTIDX = float(currline[9:])

    return STT_IMJD, STT_SMJD, PKTIDX


def calculate_start_time(imjd, smjd, pktidx):
    fs = float(3e9)     # ADC sample rate in Hz
    ts = 1 / fs         # ADC sample time in seconds
    tc = 1024. * ts     # Coarse channel sample time in seconds
    tp = 32. * tc       # Seconds per packet (32 spectra per packet)

    second_at_start_of_block = smjd + pktidx * tp

    mjd_of_block = imjd + second_at_start_of_block/86400.0

    return second_at_start_of_block, mjd_of_block


if __name__ == "__main__":
    fname = sys.argv[1]
    sec, mjd = calculate_start_mjd(*read_header(fname))

    print "  file : {0}".format(fname)
    print "start time : {0} s".format(sec)
    print "start  MJD : {0}".format(mjd)


