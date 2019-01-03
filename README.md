# Nanomechanical analysis of AFM images (Force Volume or PeakForce QNM data)

The scripts within this project analyse both force and height data rendered from various atomic force microscopy (AFM) imaging modes. They were designed to nanomechanically characterise nuclear pore complexes: nanoscale pores that reside in the nuclear membrane of eukaryotic cells. The results from this work have been published and can be seen here: http://www.life-science-alliance.org/content/1/4/e201800142.

These scripts are desined to calculate nanomechanical data at the picoNewton force regime and nanometre length scale for AFM data obtained on soft, deformable samples, in solution. They therefore may be of use to atomic force microscopists or mechanobiologists working on the nanomechanical properties of biological samples. (They should also work for hard samples, in liquid or air).

The project is separated into 5 main scripts:

1. **NM_Load_Hertz_Save.m**
2. **NM_Select_Crop.m**
3. **NM_Collate.m**
4. **NM_Rotate_Average.m**
5. **NM_Figures.m**


## Script 1: NM_Load_Hertz_Save.m ##

Loads a series of .spm image files and their concomitant .pfc files. The deflection sensitivity (nm/V) and spring constant (N/m) of the calibrated cantilever, used during the experiment, must be entered. Performs a 1st order background subtraction on the image file. It then loads all the force curves from the .pfc file, finds the contact point for each force curve (see **FC_ContactPoint_Determination.m**), and then applies the Hertz model from the contact point to the maximum indentation. Each flattened image, matrix of force curves, and matrix of effective Young's modulus (E_eff) values (from the Hertz model) are then saved out as a data structure.

Output: Flattened image, matrix of force curves aligned by baseline and contact point, matrix of calculated E_eff values.

## Script 2: NM_Select_Crop.m ##

## Script 3: NM_Collate.m ##

## Script 4: NM_Rotate_Average.m ##

## Script 5: NM_Figures.m ##
