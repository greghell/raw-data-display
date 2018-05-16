import sys
import math
import os

fread = open(sys.argv[1],'rb')  # open first file of data set
currline = str(fread.read(80))          # reads first line
print currline
nHeaderLines = 1

while currline[0:3] != 'END':           # until reaching end of header
        currline = str(fread.read(80))  # read new header line
        print currline
        if currline[0:9] == 'OBSFREQ =':        # read cenral frequency
                cenfreq = float(currline[9:])
        if currline[0:9] == 'OBSBW   =':        # read bandwidth
                obsbw = float(currline[9:])
        if currline[0:9] == 'BLOCSIZE=':        # read block size
                nblocsize = float(currline[9:])
        if currline[0:9] == 'DIRECTIO=':        # read directio flag
                ndirectio = float(currline[9:])
        nHeaderLines = nHeaderLines + 1         # counts number of lines in header

fread.close()
nPadd = 0
if ndirectio == 1:
        nPadd = int((math.floor(80.*nHeaderLines/512.)+1)*512 - 80*nHeaderLines)
statinfo = os.stat(sys.argv[1])
NumBlocks = int(round(statinfo.st_size / (nblocsize + nPadd + 80*nHeaderLines)))

FreqLow = cenfreq - obsbw/2.
FreqUp = cenfreq + obsbw/2.
FreqMin = min(FreqLow,FreqUp)
FreqMax = max(FreqLow,FreqUp)

print ""
print "minimum frequency : " + str(FreqMin) + " MHz"
print "maximum frequency : " + str(FreqMax) + " MHz"
print str(NumBlocks) + " blocks in the file"