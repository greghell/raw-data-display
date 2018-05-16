# plot_spectrogram_multi_coarse.py fname nChanOI nResol nPol

import glob
import os
import sys
import math
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) < 5:
        print '\nmissing information\n'
        print 'usage is:\n'
        print 'plot_spectrogram_coarse.py fname1 fname2 ... nChanOI nResol nPol'
        sys.exit()

filenames = []
for val in range(1,len(sys.argv)-3):
        filenames.append(sys.argv[val])

nChanOI = int(sys.argv[-3])
nResol = int(sys.argv[-2])
nPol = int(sys.argv[-1])

if (nPol != 0) and (nPol != 1):
        print 'Polarization index should be 0 or 1. Exiting...\n'
        sys.exit()

CumIntTime = np.zeros((len(filenames)))
nDataSet = 0
SpectrogTot = np.array([]).reshape((0,nResol))
for fileinit in filenames:

        print 'processing dataset #' + str(nDataSet)

        if (fileinit[-4:] == '.raw'):
                fileinit = fileinit.replace(fileinit[len(fileinit)-8:],"*.raw")
        else:
                print 'enter raw file name\n'
                sys.exit()

        all_filenames = sorted(glob.glob(fileinit))

        fname = all_filenames[0]
        fread = open(fname,'rb')        # open first file of data set
        currline = str(fread.read(80))          # reads first line
        nHeaderLines = 1
        while currline[0:3] != 'END':           # until reaching end of header
                currline = str(fread.read(80))  # read new header line
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

        for nChan in range(int(obsnchan)):
                fLowChan = fLow + (nChan)*dChanBW
                fHighChan = fLowChan + dChanBW
                print 'channel #' + str(nChan) + ' -> ' + str(min(fLowChan,fHighChan)) + ' - ' + str(max(fLowChan,fHighChan))

        NumBlockTotal = 0
        BlocksPerFile = np.zeros((len(all_filenames)))
        idx = 0
        for fname in all_filenames:
                statinfo = os.stat(fname)
                BlocksPerFile[idx] = int(statinfo.st_size / (nblocsize + nPadd + 80*nHeaderLines)/4)
                idx = idx+1

        NumBlockTotal = sum(BlocksPerFile)

        Spectrog = np.zeros((int(NumBlockTotal),int(nResol)))

        idx = 0
        nRow = 0
        for fname in all_filenames:
                print 'processing file ' + fname
                fread = open(fname,'rb')        # open first file of data set
                for nBlock in range(int(BlocksPerFile[idx])):
                        specavg = [0]*nResol
                        for nSubSpec in range(4):
                                fread.seek((4*nBlock+nSubSpec)*(nHeaderLines*80+nPadd+nblocsize)+nHeaderLines*80+nPadd+nChanOI*nChanSize)
                                sig = np.fromfile(fread, dtype = np.int8, count = nChanSize)
                                sig = np.reshape(sig,(2,nChanSize/2), order='F')
                                sigCom = np.zeros(nChanSize/2,dtype=np.complex)
                                sigCom.real = sig[0,:]
                                sigCom.imag = sig[1,:]
                                sig = np.reshape(sigCom,(2,nChanSize/4), order='F')

                                spec = np.reshape(sig[nPol,0:int(sig.shape[1]/nResol)*nResol],(nResol,int(sig.shape[1]/nResol)), order='F')
                                spec = np.power(abs(np.fft.fft(spec,axis=0)),2)

                                specavg = specavg + np.mean(spec,axis=1) / 4
                        Spectrog[nRow,:] = specavg
                        nRow = nRow + 1
                fread.close()
                idx = idx+1
        CumIntTime[nDataSet] = float(nChanSize*NumBlockTotal / (abs(obsbw)*1000000/obsnchan))
        SpectrogTot = np.concatenate((SpectrogTot,Spectrog),axis = 0)
        nDataSet = nDataSet + 1


fLow = cenfreq - obsbw/2.
fHigh = cenfreq + obsbw/2.
dChanBW = obsbw/obsnchan
fLowChan = fLow + (nChanOI)*dChanBW
fHighChan = fLowChan + dChanBW
TotInt = np.sum(CumIntTime)

plt.imshow(10*np.log10(np.flipud(SpectrogTot)),aspect='auto',extent=[fLowChan,fHighChan,0,TotInt])
plt.xlabel('frequency [MHz]',fontsize=20)
plt.ylabel('time [s]',fontsize=20)
cbar = plt.colorbar()
ax = plt.gca()
cbar.ax.tick_params(labelsize=16)
cbar.ax.set_ylabel('Power [dB]', rotation=270, fontsize=20)
plt.tick_params(labelsize=16)
for val in range(len(CumIntTime)):
        plt.plot((fLowChan,fHighChan),(np.sum(CumIntTime[0:val+1]),np.sum(CumIntTime[0:val+1])),'m')
        plt.text(fLowChan+0.1*(fHighChan-fLowChan),np.sum(CumIntTime[0:val+1])-CumIntTime[val],filenames[val],fontsize=20,color='w')
plt.title('observation time : ' + str(TotInt) + ' s',fontsize=20)
plt.grid(True)
plt.show()