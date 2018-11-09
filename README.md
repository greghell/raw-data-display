# raw-data-display  
  
# plot_coarse_channel.py  
usage:  
plot_coarse_channel.py fname fFreq nResol nPol  
plots the PSD of a coarse channel over a complete dataset around frequency fFreq, with resolution nResol for polarization nPol

# plot_coarse_channel_1file.py  
usage:  
plot_coarse_channel_1file.py fname fFreq nResol nPol  
plots the PSD of a coarse channel over a single raw file around frequency fFreq, with resolution nResol for polarization nPol

# plot_spectrogram_coarse.py
usage:  
plot_spectrogram_coarse.py fname nChanOI nResol nPol
plots the spectrogram of a given coarse channel nChanOI over a single dataset, frequency resolution nResol and polarization # nPol

# plot_coarse_spectrogram_multi.py
usage:  
plot_spectrogram_coarse.py fname0 fname1 fname2 ... nChanOI nResol nPol  
plots the spectrogram of a given coarse channel nChanOI over multiple datasets to which each of the given fname belong, frequency resolution nResol and polarization # nPol

# header_read.py
usage:  
header_read.py  
displays the raw header of file fname and gives some info on frequency span  

# plot_block.py
usage:  
plot_block.py fname nBlocOI nResol nPol  
plots spectrum of all coarse channels of a given data block nBlocOI at nResol-points resolution per coarse channel, polarization nPol

# start_time_calc.py  
usage:  
start_time_calc.py fname.raw
computes the starting time in seconds and MJD of the first block of the RAW file.
