# Nanomechanical analysis of AFM images (Force Volume or PeakForce QNM data)

These scripts analyse both force and height data rendered from various atomic force microscopy (AFM) imaging modes. They are designed to work with .spm and .pfc files produced from Bruker AFM systems. Bruker's MATLAB toolbox (NSMatlabUtilities) must therefore be installed to upload these file types into MATLAB.

The project is separated into 5 main scripts:

1. **NM_Load_Hertz_Save.m**

Applied Hertz model to each force curve in image and rotationally averaged the force data
