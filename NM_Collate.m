%%%=== NM_Collate ===%%%
 
% This is script 3/5 for the nanomechanical analysis procedure.

% This script loads all the data structures containing the cropped pores
% and all their concomitant information. It concatonates everything, and
% saves the data into a new data structure with three fields: 1, with the
% information on the cropped pores in matrices; 2, with the information on
% the cropped pores stored in their radially binned format; and 3, with the 
% information stored as the original images (uncropped). This therefore c
% concatonates all the information for both rotational averaging, and plotting.

% The resultant data structure is: 
% GenericFileName - Height_YM_CP_Cropped_rb_concatonated.mat

%% Enter load and save directories and file names etc

clear variables
close all
clc

display('NM_Collate')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%=== USER INPUT ===%%%%%%%%

% enter file numbers that would like to aggregate
FileNos = [40, 41];

%%%%%%%=== File to be loaded === %%%%%%%%%%%%%%%%%
GenericFileName = '2kHz_test_cyto';

%%%%%%%=== Save name === %%%%%%%%%%%%%%%%%
SaveName        = '2kHz_test_cyto';
%%%%%%%=== File to be loaded === %%%%%%%%%%%%%%%%%

%%%%%%%=== Data structure to be loaded
LoadFolder = ['Z:\Users\George\Documents\PhD\Data\'...
    'Nanomechanical_Outputs_Hertz_YM_CP\Test'];

%%%%%%%=== Output folder
OutputFolder = ['Z:\Users\George\Documents\PhD\Data\'...
    'Nanomechanical_Outputs_Hertz_YM_CP\Test'];
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%=== END USER INPUT ===%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create file names to be loaded

Height_YM_CP_Cropped_rb_data_structure = ' - Height_YM_CP_Cropped_rb_';
File_Suffix = strcat(GenericFileName, Height_YM_CP_Cropped_rb_data_structure);
Save_Suffix = strcat(SaveName, Height_YM_CP_Cropped_rb_data_structure);
FileType = '.mat';

[FileNames] = NanoMechFuncs.Create_FileNames_Cell_3_any(LoadFolder, File_Suffix, FileNos, FileType);


%% Store the information from the first file. Subsequent files will be concatonated to this

display('Collating data...')

load(FileNames{1});

heightdata_nm_cell_cat     = cell(1, 1);
YM_MPa_cell_cat            = cell(1, 1);
FC_nN_nm_mat_cell_cat      = cell(1, 1);
CP_Idx_Matrix_cell_cat     = cell(1, 1);
FC_count_cell_cat          = cell(1, 1);
CircleMatrix_cell_cat      = cell(1, 1);
centres_cropped_idx_cat    = cell(1, 1);
centres_cropped_xy_cat     = cell(1, 1);
centres_entire_img_idx_cat = cell(1, 1);

% store cropped matrix data
heightdata_nm_cell_cat     = Height_YM_CP_Cropped_rb.matrices.heightdata_nm_cell;
YM_MPa_cell_cat            = Height_YM_CP_Cropped_rb.matrices.YM_MPa_cell;
FC_nN_nm_mat_cell_cat      = Height_YM_CP_Cropped_rb.matrices.FC_nN_nm_mat_cell;
CP_Idx_Matrix_cell_cat     = Height_YM_CP_Cropped_rb.matrices.CP_Idx_Matrix_cell;
FC_count_cell_cat          = Height_YM_CP_Cropped_rb.matrices.FC_count_cell;
CircleMatrix_cell_cat      = Height_YM_CP_Cropped_rb.matrices.CircleMatrix_cell;
centres_cropped_idx_cat{1} = Height_YM_CP_Cropped_rb.matrices.centres_cropped_idx;
centres_cropped_xy_cat{1}  = Height_YM_CP_Cropped_rb.matrices.centres_cropped_xy;

centres_cropped_idx_cat_cell = cell(1, length(centres_cropped_idx_cat{1}(:,end)));
centres_cropped_xy_cat_cell = cell(1, length(centres_cropped_xy_cat{1}(:,end)));
for i = 1:length(centres_cropped_idx_cat{1}(:,end))
    centres_cropped_idx_cat_cell{i} = centres_cropped_idx_cat{1}(i,:);
    centres_cropped_xy_cat_cell{i}  = centres_cropped_xy_cat{1}(i,:);
end


% store radially binned data
Radial_Values_cell_cat       = Height_YM_CP_Cropped_rb.radiallybinned.Radial_Values_cell;
Circle_coords_cell_cat       = Height_YM_CP_Cropped_rb.radiallybinned.Circle_coords_cell;
heightdata_nm_rd_cropped_cat = Height_YM_CP_Cropped_rb.radiallybinned.heightdata_nm_rd_cropped;
YM_Matrix_MPa_rd_cropped_cat = Height_YM_CP_Cropped_rb.radiallybinned.YM_Matrix_MPa_rd_cropped;
FC_nN_nm_mat_rd_cropped_cat  = Height_YM_CP_Cropped_rb.radiallybinned.FC_nN_nm_mat_rd_cropped;
CP_Idx_Matrix_rd_cropped_cat = Height_YM_CP_Cropped_rb.radiallybinned.CP_Idx_Matrix_rd_cropped;
FC_count_mat_rd_cropped_cat  = Height_YM_CP_Cropped_rb.radiallybinned.FC_count_mat_rd_cropped;


% store data on entire images - first, pre-allocate where possible (this
% also forces data to be saved as cells, rather than huge, concatonated arrays).
heightdata_nm_cat = cell(1, length(FileNames));
YM_Matrix_MPa_cat = cell(1, length(FileNames));
FC_nN_nm_mat_cat  = cell(1, length(FileNames));
CP_Idx_Matrix_cat = cell(1, length(FileNames));
FC_count_mat_cat  = cell(1, length(FileNames));

heightdata_nm_cat{1}             = Height_YM_CP_Cropped_rb.entirematrices.heightdata_nm;
YM_Matrix_MPa_cat{1}             = Height_YM_CP_Cropped_rb.entirematrices.YM_Matrix_MPa;
FC_nN_nm_mat_cat{1}              = Height_YM_CP_Cropped_rb.entirematrices.FC_nN_nm_mat;
CP_Idx_Matrix_cat{1}             = Height_YM_CP_Cropped_rb.entirematrices.CP_Idx_Matrix;
FC_count_mat_cat{1}              = Height_YM_CP_Cropped_rb.entirematrices.FC_count_mat;
centres_entire_img_idx_cat{1}    = Height_YM_CP_Cropped_rb.entirematrices.centres_entire_img_idx;

%% go through all the subsequent file names and concatonate the data

for i = 1:length(FileNames)-1
    
    load(FileNames{i+1});
    % horizontally concatonate all the matrix data
    hd_cell_cat    = Height_YM_CP_Cropped_rb.matrices.heightdata_nm_cell;
    YM_cell_cat    = Height_YM_CP_Cropped_rb.matrices.YM_MPa_cell;
    FC_cell_cat    = Height_YM_CP_Cropped_rb.matrices.FC_nN_nm_mat_cell;
    CP_cell_cat    = Height_YM_CP_Cropped_rb.matrices.CP_Idx_Matrix_cell;
    count_cell_cat = Height_YM_CP_Cropped_rb.matrices.FC_count_cell;
    CM_cell_cat    = Height_YM_CP_Cropped_rb.matrices.CircleMatrix_cell;
    cc_idx_cat     = Height_YM_CP_Cropped_rb.matrices.centres_cropped_idx;
    cc_xy_cat      = Height_YM_CP_Cropped_rb.matrices.centres_cropped_xy;
    
    cc_idx_cat_cell = cell(1, length(cc_idx_cat(:,end)));
    cc_xy_cat_cell = cell(1, length(cc_xy_cat(:,end)));
    for j = 1:length(cc_idx_cat(:,end))
        cc_idx_cat_cell{j} = cc_idx_cat(j,:);
        cc_xy_cat_cell{j}  = cc_xy_cat(j,:);
    end
    
    heightdata_nm_cell_cat       = horzcat(heightdata_nm_cell_cat, hd_cell_cat);
    YM_MPa_cell_cat              = horzcat(YM_MPa_cell_cat, YM_cell_cat);
    FC_nN_nm_mat_cell_cat        = horzcat(FC_nN_nm_mat_cell_cat, FC_cell_cat);
    CP_Idx_Matrix_cell_cat       = horzcat(CP_Idx_Matrix_cell_cat, CP_cell_cat);
    FC_count_cell_cat            = horzcat(FC_count_cell_cat, count_cell_cat);
    CircleMatrix_cell_cat        = horzcat(CircleMatrix_cell_cat, CM_cell_cat);
    centres_cropped_idx_cat      = horzcat(centres_cropped_idx_cat, cc_idx_cat);
    centres_cropped_xy_cat       = horzcat(centres_cropped_xy_cat, cc_xy_cat);
    centres_cropped_idx_cat_cell = horzcat(centres_cropped_idx_cat_cell, cc_idx_cat_cell);
    centres_cropped_xy_cat_cell  = horzcat(centres_cropped_xy_cat_cell, cc_xy_cat_cell);
    
    % horizontally concatonate all the radially binned data
    RV_rb_cat     = Height_YM_CP_Cropped_rb.radiallybinned.Radial_Values_cell;
    coords_rb_cat = Height_YM_CP_Cropped_rb.radiallybinned.Circle_coords_cell;
    hd_rb_cat     = Height_YM_CP_Cropped_rb.radiallybinned.heightdata_nm_rd_cropped;
    YM_rb_cat     = Height_YM_CP_Cropped_rb.radiallybinned.YM_Matrix_MPa_rd_cropped;
    FC_rb_cat     = Height_YM_CP_Cropped_rb.radiallybinned.FC_nN_nm_mat_rd_cropped;
    CP_rb_cat     = Height_YM_CP_Cropped_rb.radiallybinned.CP_Idx_Matrix_rd_cropped;
    count_rb_cat  = Height_YM_CP_Cropped_rb.radiallybinned.FC_count_mat_rd_cropped;
    
    Radial_Values_cell_cat       = horzcat(Radial_Values_cell_cat, RV_rb_cat);
    Circle_coords_cell_cat       = horzcat(Circle_coords_cell_cat, coords_rb_cat);
    heightdata_nm_rd_cropped_cat = horzcat(heightdata_nm_rd_cropped_cat, hd_rb_cat);
    YM_Matrix_MPa_rd_cropped_cat = horzcat(YM_Matrix_MPa_rd_cropped_cat, YM_rb_cat);
    FC_nN_nm_mat_rd_cropped_cat  = horzcat(FC_nN_nm_mat_rd_cropped_cat, FC_rb_cat);
    CP_Idx_Matrix_rd_cropped_cat = horzcat(CP_Idx_Matrix_rd_cropped_cat, CP_rb_cat);
    FC_count_mat_rd_cropped_cat  = horzcat(FC_count_mat_rd_cropped_cat, count_rb_cat);
    
    % horizontally concatonate all the entire matrices
    hd_em_cat     = Height_YM_CP_Cropped_rb.entirematrices.heightdata_nm;
    YM_em_cat     = Height_YM_CP_Cropped_rb.entirematrices.YM_Matrix_MPa;
    FC_em_cat     = Height_YM_CP_Cropped_rb.entirematrices.FC_nN_nm_mat;
    CP_em_cat     = Height_YM_CP_Cropped_rb.entirematrices.CP_Idx_Matrix;
    count_em_cat  = Height_YM_CP_Cropped_rb.entirematrices.FC_count_mat;
    cc_em_idx_cat = Height_YM_CP_Cropped_rb.entirematrices.centres_entire_img_idx;
    
    cc_em_idx_cat_cell = cell(1, length(cc_em_idx_cat(:,end)));
    for j = 1:length(cc_em_idx_cat(:,end))
        cc_em_idx_cat_cell{j} = cc_em_idx_cat(j,:);
    end
    
    heightdata_nm_cat{i+1}          = Height_YM_CP_Cropped_rb.entirematrices.heightdata_nm;
    YM_Matrix_MPa_cat{i+1}          = Height_YM_CP_Cropped_rb.entirematrices.YM_Matrix_MPa;
    FC_nN_nm_mat_cat{i+1}           = Height_YM_CP_Cropped_rb.entirematrices.FC_nN_nm_mat;
    CP_Idx_Matrix_cat{i+1}          = Height_YM_CP_Cropped_rb.entirematrices.CP_Idx_Matrix;
    FC_count_mat_cat{i+1}           = Height_YM_CP_Cropped_rb.entirematrices.FC_count_mat;
    centres_entire_img_idx_cat      = horzcat(centres_entire_img_idx_cat, cc_em_idx_cat);
    
end


%%  Store in new data structure, and save

display('Saving and exporting...')

% Create a writable .mat file, and add data structures to it (to save out).
% Have to save out this way because data structures can be very large (and
% if too large cannot export using usual MATLAB save function).
FullFileOutput              = fullfile(OutputFolder, strcat(Save_Suffix, 'concatonated'));
Height_YM_CP_Cropped_rb_cat = matfile(strcat(FullFileOutput, '.mat'), 'Writable', true);

% cropped matrices
matrices.heightdata_nm_cell_cat  = heightdata_nm_cell_cat;
matrices.YM_MPa_cell_cat         = YM_MPa_cell_cat;
matrices.FC_nN_nm_mat_cell_cat   = FC_nN_nm_mat_cell_cat;
matrices.CP_Idx_Matrix_cell_cat  = CP_Idx_Matrix_cell_cat;
matrices.FC_count_cell_cat       = FC_count_cell_cat;
matrices.CircleMatrix_cell_cat   = CircleMatrix_cell_cat;
matrices.centres_cropped_idx_cat_cell = centres_cropped_idx_cat_cell;
matrices.centres_cropped_xy_cat_cell = centres_cropped_xy_cat_cell;
matrices.centres_cropped_idx_cat = centres_cropped_idx_cat;
matrices.centres_cropped_xy_cat  = centres_cropped_xy_cat;

% write to the matfile
Height_YM_CP_Cropped_rb_cat.matrices = matrices;

% entire matrices
entirematrices.heightdata_nm_cat               = heightdata_nm_cat;
entirematrices.YM_Matrix_MPa_cat               = YM_Matrix_MPa_cat;
entirematrices.FC_nN_nm_mat_cat                = FC_nN_nm_mat_cat;
entirematrices.CP_Idx_Matrix_cat               = CP_Idx_Matrix_cat;
entirematrices.FC_count_mat_cat                = FC_count_mat_cat;
entirematrices.centres_entire_img_idx_cat      = centres_entire_img_idx_cat;

% write to the matfile
Height_YM_CP_Cropped_rb_cat.entirematrices = entirematrices;

% radially binned data
radiallybinned.Radial_Values_cell_cat       = Radial_Values_cell_cat;
radiallybinned.Circle_coords_cell_cat       = Circle_coords_cell_cat;
radiallybinned.heightdata_nm_rd_cropped_cat = heightdata_nm_rd_cropped_cat;
radiallybinned.YM_Matrix_MPa_rd_cropped_cat = YM_Matrix_MPa_rd_cropped_cat;
radiallybinned.FC_nN_nm_mat_rd_cropped_cat  = FC_nN_nm_mat_rd_cropped_cat;
radiallybinned.CP_Idx_Matrix_rd_cropped_cat = CP_Idx_Matrix_rd_cropped_cat;
radiallybinned.FC_count_mat_rd_cropped_cat  = FC_count_mat_rd_cropped_cat;

% write to the matfile
Height_YM_CP_Cropped_rb_cat.radiallybinned = radiallybinned;


