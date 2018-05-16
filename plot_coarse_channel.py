# plot_coarse_channel.py fname fFreq nResol nPol

import glob
import os
import sys
import math
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) != 5:
        print '\nmissing information\n'
        print 'usage is:\n'
        print '''plot_coarse_channel.py fname fFreq nResol nPol'''
        sys.exit()

fileinit = str(sys.argv[1])
fFreq = float(sys.argv[2])
nResol = int(sys.argv[3])
nPol = int(sys.argv[4])

if (nPol != 0) and (nPol != 1):
        print 'Polarization index should be 0 or 1. Exiting...\n'
        sys.exit()

if (fileinit[-4:] == '.raw'):
        fileinit = fileinit.replace(fileinit[len(fileinit)-8:],"*.raw")
else:
        print 'enter raw file name\n'
        sys.exit()

all_filenames = sorted(glob.glob(fileinit))
print all_filenames

fname = all_filenames[0]
fread = open(fname,'rb')        # open first file of data set
currline = str(fread.read(80))          # reads first line
nHeaderLines = 1
while currline[0:3] != 'END':           # until reaching end of header
        currline = str(fread.read(80))  # read new header line
        # print currline
        if currline[0:9] == 'OBSFREQ =':        # read cenral frequency
                cenfreq = float(currline[9:])
        if currline[0:9] == 'OBSBW   =':        # read bandwidth
                obsbw = float(currline[9:])
        if currline[0:9] == 'OBSNCHAN=':        # read number of coarse channels
                obsnchan = float(currline[9:])
        if currline[0:9] == 'DIRECTIO=':        # read directio flag
                ndirectio = float(currline[9:])
        if currline[0:9] == 'BLOCSIZE=':        # read block size
                nblocsize = float(currline[9:])

        nHeaderLines = nHeaderLines + 1         # counts number of lines in header

fread.close()

nChanSize = int(nblocsize / obsnchan)
nPadd = 0
if ndirectio == 1:
        nPadd = int((math.floor(80.*nHeaderLines/512.)+1)*512 - 80*nHeaderLines)
statinfo = os.stat(fname)
NumBlocs = int(round(statinfo.st_size / (nblocsize + nPadd + 80*nHeaderLines)))

fLow = cenfreq - obsbw/2.
fHigh = cenfreq + obsbw/2.
dChanBW = obsbw/obsnchan

if fFreq < min(fLow,fHigh) or fFreq > max(fLow,fHigh):
		print 'Frequency not covered by file\n'
		print 'Frequency bandwidth = ['+str(min(fLow,fHigh))+','+str(max(fLow,fHigh))+']\n'
		sys.exit()

#for nChan in range(int(obsnchan)):
#        fLowChan = fLow + (nChan)*dChanBW
#        fHighChan = fLowChan + dChanBW
#        print 'channel #' + str(nChan) + ' -> ' + str(min(fLowChan,fHighChan)) + ' - ' + str(max(fLowChan,fHighChan))

if obsbw > 0:
		nChanOI = int((fFreq-fLow)/dChanBW)
else:
		nChanOI = int((fLow-fFreq)/abs(dChanBW))

specavg = [0]*nResol
NumBlockTotal = 0

for fname in all_filenames:

        print 'processing file ' + fname

# fname = sys.argv[1]
        fread = open(fname,'rb')        # open first file of data set
        currline = str(fread.read(80))          # reads first line
        nHeaderLines = 1
        while currline[0:3] != 'END':           # until reaching end of header
                currline = str(fread.read(80))  # read new header line
                # print currline
                if currline[0:9] == 'OBSFREQ =':        # read cenral frequency
                        cenfreq = float(currline[9:])
                if currline[0:9] == 'OBSBW   =':        # read bandwidth
                        obsbw = float(currline[9:])
                if currline[0:9] == 'OBSNCHAN=':        # read number of coarse channels
                        obsnchan = float(currline[9:])
                if currline[0:9] == 'DIRECTIO=':        # read directio flag
                        ndirectio = float(currline[9:])
                if currline[0:9] == 'BLOCSIZE=':        # read block size
                        nblocsize = float(currline[9:])

                nHeaderLines = nHeaderLines + 1         # counts number of lines in header

        nChanSize = int(nblocsize / obsnchan)
        nPadd = 0
        if ndirectio == 1:
                nPadd = int((math.floor(80.*nHeaderLines/512.)+1)*512 - 80*nHeaderLines)
        statinfo = os.stat(fname)
        NumBlocs = int(round(statinfo.st_size / (nblocsize + nPadd + 80*nHeaderLines)))
        NumBlockTotal = NumBlockTotal + NumBlocs

        for nBlock in range(NumBlocs):
                fread.seek(nBlock*(nHeaderLines*80+nPadd+nblocsize)+nHeaderLines*80+nPadd+nChanOI*nChanSize)
                sig = np.fromfile(fread, dtype = np.int8, count = nChanSize)
                sig = np.reshape(sig,(2,nChanSize/2), order='F')
                sigCom = np.zeros(nChanSize/2,dtype=np.complex)
                sigCom.real = sig[0,:]
                sigCom.imag = sig[1,:]
                sig = np.reshape(sigCom,(2,nChanSize/4), order='F')

                spec = np.reshape(sig[nPol,0:int(math.floor(sig.shape[1]/nResol)*nResol)],(nResol,int(math.floor(sig.shape[1]/nResol))), order='F')
                spec = np.power(abs(np.fft.fft(spec,axis=0)),2)

                specavg = specavg + np.mean(spec,axis=1) / NumBlocs

        fread.close()


fLow = cenfreq - obsbw/2.
fHigh = cenfreq + obsbw/2.
dChanBW = obsbw/obsnchan
fLowChan = fLow + (nChanOI)*dChanBW
fHighChan = fLowChan + dChanBW
TotInt = float(nChanSize/4*NumBlockTotal / (abs(obsbw)*1000000/obsnchan))

plt.plot(np.linspace(fLowChan,fHighChan,nResol), 10*np.log10(np.fft.fftshift(specavg)),'k')
plt.xlabel('frequency [MHz]',fontsize=20)
plt.ylabel('power spectral density [dB]',fontsize=20)
plt.tick_params(labelsize=16)
plt.suptitle(sys.argv[1])
plt.title('integration time : ' + str(TotInt) + ' s')
plt.grid(True)
plt.show()