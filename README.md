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

Loads the image, E_eff, and aligned force curve data previously saved into a .mat structure. Asks the user to select the centres of the pores, crops the pores and E_eff values, and saves the cropped image and E_eff data into new matrices. It then radially bins this data.

This script is designed to work for pores cropped near the edge of the image. As long as the central axis of rotation is visible to the user, it can be selected, and the data will be carried forward.

The output is a new data structure containing the cropped height and E_eff data, along with the radially binned height and E_eff data. A count array comes with each radially binned array, tracking how many force curves are in each radial bin. This takes into account bins of different sizes, and force curves that have been binned.

## Script 3: NM_Collate.m ##

This script is designed to aggregate all the analysed images within one experiment and put them into one data structure. This enables averaging of the results in the following script.

Loads all the data structures containing the cropped pores and all their concomitant information (must enter file numbers would like to concatonate, manually). It concatonates everything, and saves the data into a new data structure with three fields: 1, with the information on the cropped pores in matrices; 2, with the information on the cropped pores stored in their radially binned format; and 3, with the information stored as the original images (uncropped). This therefore concatonates all the information for both rotational averaging, and plotting.

## Script 4: NM_Rotate_Average.m ##

This script loads all the concatonated data from the previous script. It then averages all the data based on their radial distance from the centre of the pore - ready for plotting in the next (and final) script.

## Script 5: NM_Figures.m ##

This script generates publication quality figures.

## Authors

George J Stanley

## Acknowledgements

Acknowledgements to Jonathan Lanset for linspecer, and to Rob Campbell for shadedErrorBar.
