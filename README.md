# android-denoise-algorithm
Android noise-reduction algorithms.  
This class provides some noise-reduction algorithms based on this paper: http://www.citi.sinica.edu.tw/papers/yu.tsao/4900-F.pdf  
It utilizes the Fast Fourier transform (FFT) functions from libgdx audio library, and manipulates input audio to output noise-reduced data.

## Overview
The class constructor takes the desired FFT length and sample frequency values.  
The method getSigma() takes audio data and calculates the average noise parameter (as a constant). The method denoise() takes input audio data and algorithm type and computes the output data which is noise-reduced.
