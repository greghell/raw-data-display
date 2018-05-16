# pt_block.py fname nBlocOI nResol nPol

import glob
import os
import sys
import math
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) != 5:
        print '\nmissing information\n'
        print 'usage is:\n'
        print '''plot_coarse_channel.py fname nBlocOI nResol nPol'''
        sys.exit()

fname = str(sys.argv[1])
nBlocOI = int(sys.argv[2])
nResol = int(sys.argv[3])
nPol = int(sys.argv[4])

if (nPol != 0) and (nPol != 1):
        print 'Polarization index should be 0 or 1. Exiting...\n'
        sys.exit()

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
                obsnchan = int(currline[9:])
        if currline[0:9] == 'DIRECTIO=':        # read directio flag
                ndirectio = int(currline[9:])
        if currline[0:9] == 'BLOCSIZE=':        # read block size
                nblocsize = int(currline[9:])

        nHeaderLines = nHeaderLines + 1         # counts number of lines in header

nChanSize = int(nblocsize / obsnchan)
nPadd = 0
if ndirectio == 1:
        nPadd = int((math.floor(80.*nHeaderLines/512.)+1)*512 - 80*nHeaderLines)
statinfo = os.stat(fname)
NumBlocs = int(round(statinfo.st_size / (nblocsize + nPadd + 80*nHeaderLines)))

fLow = cenfreq - obsbw/2.
fHigh = cenfreq + obsbw/2.
dChanBW = obsbw/obsnchan

for nChan in range(int(obsnchan)):
        fLowChan = fLow + (nChan)*dChanBW
        fHighChan = fLowChan + dChanBW

specavg = [0]*nResol*obsnchan
NumBlockTotal = 0

# plt.hold(True)
for nChan in range(int(obsnchan)):
        fread.seek(nBlocOI*(nHeaderLines*80+nPadd+nblocsize)+nHeaderLines*80+nPadd+nChan*nChanSize)
        sig = np.fromfile(fread, dtype = np.int8, count = nChanSize)
        sig = np.reshape(sig,(2,nChanSize/2), order='F')
        sigCom = np.zeros(nChanSize/2,dtype=np.complex)
        sigCom.real = sig[0,:]
        sigCom.imag = sig[1,:]
        sig = np.reshape(sigCom,(2,nChanSize/4), order='F')

        spec = np.reshape(sig[nPol,0:int(math.floor(sig.shape[1]/nResol)*nResol)],(nResol,int(math.floor(sig.shape[1]/nResol))), order='F')
        spec = np.power(abs(np.fft.fft(spec,axis=0)),2)

        specavg = np.mean(spec,axis=1)
        plt.plot(np.linspace(fLow + (nChan)*dChanBW,fLow + (nChan)*dChanBW + dChanBW,nResol), 10*np.log10(np.fft.fftshift(specavg)),'k')

plt.xlabel('frequency [MHz]',fontsize=20)
plt.ylabel('power spectral density [dB]',fontsize=20)
plt.tick_params(labelsize=16)
plt.grid(True)
plt.show()
# plt.show()