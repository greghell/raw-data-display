# calculates the start time in seconds and MJD of a given RAW file.
# call : python start_time_calc.py filename.raw
# (code follows https://github.com/UCBerkeleySETI/breakthrough/blob/master/doc/Time-in-RAW-Files.md)

import sys
import math
import os

fread = open(sys.argv[1],'rb')                  # open first file of data set
currline = str(fread.read(80))          # reads first line

while currline[0:3] != 'END':           # until reaching end of header
        currline = str(fread.read(80))  # read new header line
        if currline[0:9] == 'STT_IMJD=':        # read cenral frequency
                STT_IMJD = float(currline[9:])
        if currline[0:9] == 'STT_SMJD=':        # read bandwidth
                STT_SMJD = float(currline[9:])
        if currline[0:9] == 'PKTIDX  =':        # read block size
                PKTIDX = float(currline[9:])

fread.close()

fs = float(3e9) # ADC sample rate in Hz
ts = 1/fs       # ADC sample time in seconds
tc = 1024. * ts # Coarse channel sample time in seconds
tp = 32.*tc             # Seconds per packet (32 spectra per packet)

second_at_start_of_block = STT_SMJD + PKTIDX * tp

mjd_of_block = STT_IMJD + second_at_start_of_block/86400.0

print ""
print "start time : " + str(second_at_start_of_block) + " s"
print "start MJD : " + str(mjd_of_block)
