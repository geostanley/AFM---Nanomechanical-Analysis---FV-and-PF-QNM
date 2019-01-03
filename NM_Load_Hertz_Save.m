%%%=== NM_Load_Hertz_Save ===%%%

% This is script 1/5 for the nanomechanical analysis procedure.

% It loads a series of .spm image files and their concomitant .pfc files.
% It performs a 1st order background subtraction on the image file. It then
% loads all the force curves from the .pfc file, finds the contact point
% for each force curve and then applied the Hertz model from the contact
% point to the maximum indentation. Each flattened image, matric of force
% curves, and matrix of E_{eff} values (from the Hertz model) are then
% saved as a data structure.


%% Enter filename and directory for data

clear variables
close all
clc

display('NM_Load_Hertz_Save')

%%%=== Must fill in these parameters ===%%%

% 09/01/2018 - Sample 3 - A027-01-05 - 2kHz

% Input folder with data
Data_Directory = ['Z:\Users\George\Documents\PhD\Data\AFM\ICON\'...
    '2018_01_09 - Icon - NIM - MSNL E - PFQNM\Sample 3\A027-01-05 - MSNL E\2kHz\'];

GenericFileName = 'NIM - MSNL E - Cyt - 45mV';

% Output directory and filename
OutputFolder = ['Z:\Users\George\Documents\PhD\Data\'...
    'Nanomechanical_Outputs_Hertz_YM_CP\Test'];

Save_name = '2kHz_test_cyto';

%%%=== Enter calibration details of cantilever ===%%%
DeflSens    = 55.30; % in nm/V
SpringConst = 0.1594; % in N/m

channel = 1; % from which channel do you want to take the height data (trace or retrace, height or height sensor)?

HertzIndentation_nm = 20; % distance (nm) over which fit the Hertz model
NS_nine_one = 'no'; % if files from Nanoscope 9.1: 'yes'; if 9.2 or later: 'no'


% The correctly numbered spm file does not
% always save with the correctly numbered pfc file. Therefore, the user must
% enter the numbers of the spm and pfc files manually, ensuring that the
% correct pairs are at the same index position in the arrays. If the two
% arrays are not of the same length, an error will appear. If the pairing
% of numbers is incorrect, the Hertz model will be applied to force curves
% from the wrong image.

% Later, the saved files are numbered 'correctly', always using the spm
% file number. This is because there are usually more spm files than pfc
% files (because some pfc files do not save properly... because of bugs in 
% Bruker's software).
File_Nos_spm = [40,41];
File_Nos_pfc = [39,40];

% if PF QNM data: 1; if FV: 0
QNM = 1;


%% Load data, flatten height images, organise force curves, find contact point, and apply Hertz model

% Open NSMatlabUtilities
NSMU = NSMatlabUtilities();

if length(File_Nos_spm) ~= length(File_Nos_pfc);
    error('Number of spm files does not match number of pfc files.')
end

% Create strings of the files to be loaded
[FullFile_spm, FullFile_pfc] = NanoMechFuncs.Create_FileNames_Cell_3_spm_pfc(Data_Directory, GenericFileName, File_Nos_spm, File_Nos_pfc, NS_nine_one);

for n = 1:length(FullFile_spm)

    FullFileSPM = FullFile_spm{n};
    FullFilePFC = FullFile_pfc{n};
    
    if QNM == 0;
        FullFilePFC = FullFileSPM;
    end

    % Load the height data
    NSMU.Open(FullFileSPM);

    % Ensure the correct channel has been selected (see beginning of script)
    if QNM == 1;
        [height_mat, scale_units, type_desc] = NSMU.GetImageData(channel, NSMU.METRIC);
    else
        [height_mat, scale_units, type_desc] = NSMU.GetForceVolumeImageData(NSMU.METRIC);
    end

    % display name of channel and format of data loaded (i.e. the unit)
    disp(['Channel = ', type_desc])
    disp(['Data Unit = ', scale_units])

    [ScanSize_nm, ScanSizeUnit] = NSMU.GetScanSize(channel);
    disp(['Image Size = ',num2str(ScanSize_nm), ' ',ScanSizeUnit])


    %% Flatten the height data (1st order plane background subtraction)
    
    % transform to XYZ array for plane fitting
    [XYZ_array_1] = ImageFlatteningFuncs.Matrix_to_Nx3array(height_mat);

    % find 1sr order plane of entire image
    BinWidth_nm = 0.1;
    Plane_fit_mask_1 = 1;
    greater_than_1   = 0;
    % index XYZ array
    [XYZ_array_for_plane_fit_1, ~] = ImageFlatteningFuncs.XYZarray_indexed_by_percentage_height(XYZ_array_1, BinWidth_nm, Plane_fit_mask_1, greater_than_1);
    % get plane and subtract from height data
    [plane_1] = ImageFlatteningFuncs.PlaneFit_XYZarray(height_mat, XYZ_array_for_plane_fit_1);
    height_mat_2 = height_mat - plane_1;


    % now repeat but only apply plane fitting to bottom 50% of data (background)
    [XYZ_array_2] = ImageFlatteningFuncs.Matrix_to_Nx3array(height_mat_2);

    Plane_fit_mask_1 = 0.5;
    greater_than_1   = 0;

    [XYZ_array_for_plane_fit_2, ~] = ImageFlatteningFuncs.XYZarray_indexed_by_percentage_height(XYZ_array_2, BinWidth_nm, Plane_fit_mask_1, greater_than_1);
    [plane_2] = ImageFlatteningFuncs.PlaneFit_XYZarray(height_mat_2, XYZ_array_for_plane_fit_2);
    heightdata_nm = height_mat_2 - plane_2;

    %% Load force curves and save in matrix with same coordinate system as image data

    % close the spm file
    NSMU.Close();
    % open the pfc file
    NSMU.Open(FullFilePFC);
    NumberOfCurves = NSMU.GetNumberOfForceCurves();
    NumberOfLines  = NSMU.GetSamplesPerLine(1);

    % Need to make a curve number matrix to load the force curves. This needs
    % to be in the same coordinate system as the heightdata. The CurveMatrix is
    % simply a matrix of intergers monotonically increasing from 1 to the number
    % of force curves in the image. However, these force curves are not in the
    % same coordinate system as the height data: they are effectively rotated
    % 90 degrees clockwise. Therefore, to correct for this, the CurveMatrix is
    % transposed and then flipped upside down. If you run the Hertz model on an
    % entire forcurve matrix, after doing this operation, you will see the
    % resultant YoungsModuls_MPa_Matrix matches up with the height data.

    Curve_Number_Array = 1:NumberOfCurves;
    CurveMatrix_wrong_orientation = reshape(Curve_Number_Array, NumberOfLines, NumberOfLines);
    % CurveMatrix is now in the correct orientation
    CurveMatrix = flipud(CurveMatrix_wrong_orientation');

    % This section creates a Curve Matrix (just a matrix of intergers
    % monotonically increasing from 1 to the number of force curves), then uses
    % these 
    FC_Matrix = cell(size(CurveMatrix));

    [row, col] = size(CurveMatrix);

    % display name of channel and format of data loaded (i.e. the unit)
    disp('Loading force curves and storing in matrix...')

    for i = 1:row

        for j = 1:col

            CurveNumber = CurveMatrix(i,j);

            if QNM == 1;
                % Obtain the Piezo displacement in nm
                [xTrace_Metric, ~, ~, ~, xLabel, ~] = NSMU.CreatePeakForceForceCurveZplot(CurveNumber, NSMU.METRIC, 2); % 2 = Z-displacement rather than tip-sample separation
                % Deflection error in V
                [~, ~, yTrace_Volts, ~, ~, yLabel]  = NSMU.CreatePeakForceForceCurveZplot(CurveNumber, NSMU.VOLTS, 2);
            else
                [xTrace_Metric, ~, ~, ~, ~, ~] = NSMU.CreateForceVolumeForceCurveZplot(CurveNumber, NSMU.METRIC, 2); % 2 = Z-displacement rather than tip-sample separation
                % Deflection error in V
                [~, ~, yTrace_Volts, ~, ~, ~]  = NSMU.CreateForceVolumeForceCurveZplot(CurveNumber, NSMU.VOLTS, 2);
            end            

            xTrace_Piezo = flipud(xTrace_Metric);
            yTrace_Defl  = yTrace_Volts.*DeflSens; % Defl Error in nm (converted from V), for tip-sample separation
            xTrace_Sep   = xTrace_Piezo + yTrace_Defl; % Z-Piezo displacement to tip-sample separation
            yTrace_nN    = yTrace_Volts.*(DeflSens*SpringConst); % deflection in Volts to force on CL

            % save force curve
            FC_nN_nm       = [xTrace_Sep,yTrace_nN];
            FC_Matrix{i,j} = FC_nN_nm;

        end

    end
    
    %% Find contact point and apply Hertz model
    
    disp('Assigning the contact point and applying the Hertz model...')

        CP_Idx_Matrix = zeros(size(CurveMatrix));
        YM_Matrix_MPa = zeros(size(CurveMatrix));
        FC_count_mat  = zeros(size(YM_Matrix_MPa));
        FC_nN_nm_mat  = cell(size(YM_Matrix_MPa));

    tic

    parfor i = 1:row
              
        display(['CP and Hertz. Force curves ', num2str(((i*row) - (row-1))), ':', num2str((i*row)), '/', num2str(NumberOfCurves)])
        
        for j = 1:col

            FC_nN_nm = FC_Matrix{i,j};   
            Factor = 3;
            [~, FC_nN_nm_Corrected, Contact_Point_Index, Indentation] = FC_ContactPoint_Determination(FC_nN_nm, Factor, 'no', 21);

            if isnan(Contact_Point_Index) == 1 % if contact point could not be found           
                YoungsMod_MPa = 0; % using the counter, the 0 can be discounted in averaging later (by adding the 0 but dividing by the sum of the counter)
                FC_count_mat(i,j) = 0; % give the counter a 0
                FC_nN_nm_mat{i,j} = nan; % put a not-a-number where the force curve should have been
                
            else % if contact point found, do Hertz, save corrected force curve, and give counter a 1   
                
                YoungsMod_MPa = HertzModel(FC_nN_nm_Corrected, Contact_Point_Index, HertzIndentation_nm);
                FC_count_mat(i,j)  = 1; % give counter a 1
                FC_nN_nm_mat{i,j}  = FC_nN_nm_Corrected; % save corrected force curve
            end

            CP_Idx_Matrix(i,j) = Contact_Point_Index;
            YM_Matrix_MPa(i,j) = YoungsMod_MPa;
            

        end

    end

    time_CP_Hertz = toc;
    display(time_CP_Hertz)

    %% Count how many force curves have been binned

    FC_count_sum = sum(FC_count_mat(:));
    FC_bin_percentage = 100 - ((FC_count_sum/NumberOfCurves)*100);
    
    %% Display the two main outputs
    
    figure(), imagesc(heightdata_nm)
    colorbar
    g = colorbar;
    ylabel(g,'Height (nm)', 'FontSize', 13)
    set(gca,'ydir','normal');
    set(gca, 'FontSize', 13)
    title('Height data (nm)')

    figure(), imagesc(YM_Matrix_MPa)
    caxis([0 4])
    colorbar
    g = colorbar;
    ylabel(g,'YM (MPa)', 'FontSize', 13)
    set(gca,'ydir','normal');
    set(gca, 'FontSize', 13)

    %% Save the flattened height data, the contact point indices, and the Young's moduli

    FileNumb = File_Nos_spm(n);
    FullFileOutput = fullfile(OutputFolder, strcat(Save_name, ' - ', ' Height_YM_CP', '_', num2str(FileNumb)));

    Height_YM_CP.heightdata_nm                     = heightdata_nm;
    Height_YM_CP.YM_Matrix_MPa                     = YM_Matrix_MPa;
    Height_YM_CP.FC_nN_nm_mat                      = FC_nN_nm_mat;
    Height_YM_CP.CP_Idx_Matrix                     = CP_Idx_Matrix;
    Height_YM_CP.FC_count_mat                      = FC_count_mat;
    Height_YM_CP.ImageParameters.ScanSize_nm       = ScanSize_nm;
    Height_YM_CP.ImageParameters.FC_bin_percentage = FC_bin_percentage;

    save(strcat(FullFileOutput, '.mat'), 'Height_YM_CP');

    % close the pfc file before staring the next loop
    NSMU.Close();
    
end
    
    







