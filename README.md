Original release date : 2/28/18

Reference #1          : "prDeep: Robust Phase Retrieval with Flexible Deep Neural Networks"

Authors               : Christopher A. Metzler, Philip Schniter, Ashok Veeraraghavan, Richard G. Baraniuk


Questions/suggestions/comments: chris.metzler@rice.edu

Primary Contents
---------------------------------------------------------------------------
Scripts:
  * PR_demo.m: Demonstrates phase retrieval from noisy coded diffraction and fourier measurements using prDeep.

Functions:
  * prDeep.m: Implementation of prDeep
  * HIO.m: Basic implementation of HIO algorithm
    
Auxiliary functions:
  * disambig2Drfft.m: Phaseless 2D fourier transforms have phase and translation ambiguities. This function accounts for them.


Packages
---------------------------------------------------------------------------
This download includes trained DnCNN denoiser weights for various noise levels.

Dependencies
---------------------------------------------------------------------------
D-AMP Toolbox (https://github.com/ricedsp/D-AMP_Toolbox) must be on your path.
MatConvNet (http://www.vlfeat.org/matconvnet/) must be compiled and on your path.
FASTA (https://github.com/tomgoldstein/fasta-matlab/) must be on your path.

Installation
---------------------------------------------------------------------------
MatConvNet must be compiled before it can be used.
