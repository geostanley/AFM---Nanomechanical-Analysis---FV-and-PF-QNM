%%%=== NM_Select_Crop ===%%%

% This is script 2/5 for the nanomechanical analysis procedure.

% It loads the Image, YM (i.e., E_eff), and CP data previously saved into a .mat
% structure. It asks the user to select the centres of the pores, crops the
% pores and Young's moduli, and saves the cropped image data and Hertz data into new matrices.
% It then radially bins this data.

% This script is designed to work for pores cropped near the edge of the
% image. As long as the central axis of rotation is visible to the user, it
% can be selected, and the data will be carried forward.

% The output is a new data structure containing the cropped height data and
% Hertz data, along with the radially binned height and Hertz data. A count
% array comes with each Radially binned array, telling how many force
% curves are in each radial bin. This takes into account bins of different
% sizes, and force curves that have been binned.

% This script is designed to be used one file at a time, so that if the
% user makes a mistake inputting the centres, it can be easily redone.

% NB: I tested whether the indexing of the force curves was correct 
% (i.e., whether they matched up with the imcropping of the other matrices) 
% by applying the Hertz model to cropped force curves, and comparing to
% the corresponding cropped_YM matrix, and they were indeed the same.
% This uses the FC_nN_nm_mat_cell (containing indexed, or cropped,
% force curve matrices), and the corresponding CP_Idx_Matrix_cell.
% I.e., everything matches up.   

%% Enter load and save directories and file names etc

clear variables
close all
clc

display('NM_Select_Crop')

%%%%%%%=== user input file number ===%%%%%%%
FileNumb = 40;


height_min_nm = -10; 
height_max_nm = 50;

%%%%%%%=== File to be loaded === %%%%%%%%%%%%%%%%%
GenericFileName = '2kHz_test_cyto';

%%%%%%%=== Data structure to be loaded
LoadFolder = ['Z:\Users\George\Documents\PhD\Data\'...
    'Nanomechanical_Outputs_Hertz_YM_CP\Test'];

%%%%%%%=== Output folder
OutputFolder = ['Z:\Users\George\Documents\PhD\Data\'...
    'Nanomechanical_Outputs_Hertz_YM_CP\Test'];

%%%=== Enter crop size and radial bin sizes ===%%%
CropSize_nm       = 200; % usually go 200 nm for NPCs
Initial_Radius_nm = 5; % radius of first bin
Radius_nm         = 4; % further bins in increaments of how many nm
NumberOfBins      = 25; % number of bins in total

%% Load data structure

FullFileLoad = fullfile(LoadFolder, strcat(GenericFileName, ' - ', ' Height_YM_CP', '_', num2str(FileNumb), '.mat'));

load(FullFileLoad)

heightdata_nm = Height_YM_CP.heightdata_nm;
YM_Matrix_MPa = Height_YM_CP.YM_Matrix_MPa;
FC_nN_nm_mat  = Height_YM_CP.FC_nN_nm_mat;
CP_Idx_Matrix = Height_YM_CP.CP_Idx_Matrix;
FC_count_mat  = Height_YM_CP.FC_count_mat;
ScanSize_nm   = Height_YM_CP.ImageParameters.ScanSize_nm;

%% Select pore centres and crop

% calculate size of imcrop and ensure even number of pixels 
[entireImg_row, entireImg_col] = size(heightdata_nm);
nmperpixel                     = ScanSize_nm/entireImg_row;
CropSize_pix                   = round(CropSize_nm/nmperpixel);

% imcrop crops from the coordinates c1 and c2, a distance of CropSize_pix.
% This excludes the corner pixel. Therefore, the resultant crop is of the
% size CropSize_pix+1.
% This means when deciding on CropSize_pix, need to ensure it is an even
% number, as when cropping later, this will make it odd (as it will add a
% pixel). We want an odd number of pixels to ensure the rotational
% averaging is as easy as possible.
if mod(CropSize_pix, 2) ~= 0 % if numb of pixels is odd
    CropSize_pix = CropSize_pix + 1; % make even
end

% CropSize_pix will be an even number. So divide by two and subtract from
% centres. These are the starting coordinates for the image crop.
CropSize_pix_half = floor(CropSize_pix/2);

%% Now know size to be cropped, assign centres

centres_entire_img_xy = NanoMechFuncs.InputCentres_Colourmap(heightdata_nm, 'parula');

% Pre-allocate
heightdata_nm_cell     = cell(1, length(centres_entire_img_xy(:,end)));
YM_MPa_cell            = cell(1, length(centres_entire_img_xy(:,end)));
FC_nN_nm_mat_cell      = cell(1, length(centres_entire_img_xy(:,end)));
CP_Idx_Matrix_cell     = cell(1, length(centres_entire_img_xy(:,end)));
FC_count_cell          = cell(1, length(centres_entire_img_xy(:,end)));
CircleMatrix_cell      = cell(1, length(centres_entire_img_xy(:,end)));
Radial_Values_cell     = cell(1, length(centres_entire_img_xy(:,end)));
Circle_coords_cell     = cell(1, length(centres_entire_img_xy(:,end)));
centres_entire_img_idx = zeros(size(centres_entire_img_xy));
centres_cropped_idx    = zeros(size(centres_entire_img_xy));
centres_cropped_xy     = zeros(size(centres_entire_img_xy));

% Pre-allocate for radially binned data (rd) - these will become cells of
% cells
heightdata_nm_rd_cropped = cell(1, length(centres_entire_img_xy(:,end)));
YM_Matrix_MPa_rd_cropped = cell(1, length(centres_entire_img_xy(:,end)));
FC_nN_nm_mat_rd_cropped  = cell(1, length(centres_entire_img_xy(:,end)));
CP_Idx_Matrix_rd_cropped = cell(1, length(centres_entire_img_xy(:,end)));
FC_count_mat_rd_cropped  = cell(1, length(centres_entire_img_xy(:,end)));

FC_nN_nm_radially_binned_crops = cell(1, NumberOfBins);

for i = 1:length(centres_entire_img_xy(:,end))
        
    centre_x = centres_entire_img_xy(i,1);
    centre_y = centres_entire_img_xy(i,2);
    
    % save as idx
    centre_row = centre_y;
    centre_col = centre_x;
    centres_entire_img_idx(i,1) = centre_row;
    centres_entire_img_idx(i,2) = centre_col;

    % c1 and c2 are corners in xy system for imcrop. Also below is
    % coordinates system in idx for indexing force curves
    c1 = centre_x-CropSize_pix_half;
    c2 = centre_y-CropSize_pix_half;
    
    % get these in index values for cropping FC matrix
    c1_col_idx = centre_col-CropSize_pix_half;
    c3_col_idx = c1_col_idx+CropSize_pix;
    c2_row_idx = centre_row-CropSize_pix_half;
    c4_row_idx = c2_row_idx+CropSize_pix;
    
    % if one of the index values is going to exceed the matrix dimensions,
    % change it so that it is instead at the edge of the image (this is not
    % a problem for imcrop, which is why it was originally used, however, 
    % imcrop cannot be used to index the force curves in a cell array, and
    % so regular indexing must be used). This should all be done using
    % indexing really, but it works, so...
    if c2_row_idx<1
        c2_row_idx = 1;
    end
    
    if c1_col_idx<1
        c1_col_idx = 1;
    end
    
    if c4_row_idx>entireImg_row
        c4_row_idx = entireImg_row;
    end
    
    if c3_col_idx>entireImg_col
        c3_col_idx = entireImg_col;
    end

    % imcrop crops from the coordinates c1 and c2, a distance of CropSize_pix.
    % This excludes the corner pixel. Therefore, the resultant crop is of the
    % size CropSize_pix+1.
    cropped_img           = imcrop(heightdata_nm, [c1 c2 CropSize_pix CropSize_pix]);
    cropped_YM            = imcrop(YM_Matrix_MPa, [c1 c2 CropSize_pix CropSize_pix]);
    cropped_CP_Idx_Matrix = imcrop(CP_Idx_Matrix, [c1 c2 CropSize_pix CropSize_pix]);
    cropped_FC_count_mat  = imcrop(FC_count_mat, [c1 c2 CropSize_pix CropSize_pix]);
    cropped_FC_nN_nm      = FC_nN_nm_mat(c2_row_idx:c4_row_idx, c1_col_idx:c3_col_idx);
    
    %% with the cropped image, planefit the top of the pore and bring to
    % 40nm (such that all pores are at the same height for rotational
    % averaging).
    
    % transform to XYZ array for plane fitting
    [XYZ_array_1] = ImageFlatteningFuncs.Matrix_to_Nx3array(cropped_img);

    % for plane fitting, take the top 50% of the data
    BinWidth_nm = 0.1; % the smaller the value the more accurate the XYZ indexing. 0.1nm should be fine.
    Plane_fit_mask_1 = 0.50;
    greater_than_1   = 1;
    % create new XYZ array of only bottom 50% of data
    [XYZ_array_for_plane_fit_1, ~] = ImageFlatteningFuncs.XYZarray_indexed_by_percentage_height(XYZ_array_1, BinWidth_nm, Plane_fit_mask_1, greater_than_1);

    % find the plane of the indexed bottom 50% of height data
    [plane_1] = ImageFlatteningFuncs.PlaneFit_XYZarray(cropped_img, XYZ_array_for_plane_fit_1);
    % subtract plane from data. This operation also brings the rim of the
    % pore to height_rim_pore_nm (usually set to 35 nm).
    height_rim_pore_nm = 35;
    heightdata_cropped_nm = (cropped_img - plane_1) + height_rim_pore_nm;
    
    %%
    
    % NB: I tested whether the indexing of the force curves was correct by
    % applying the Hertz model to cropped force curves, and comparing to
    % the corresponding cropped_YM matrix, and they were indeed the same.
    % This uses the FC_nN_nm_mat_cell (containing indexed, or cropped,
    % force curve matrices), and the corresponding CP_Idx_Matrix_cell.
    % I.e., everything matches up.    
    
    % Save cropped height and YM data
    heightdata_nm_cell{i} = heightdata_cropped_nm;
    YM_MPa_cell{i}        = cropped_YM;
    FC_count_cell{i}      = cropped_FC_count_mat;
    CP_Idx_Matrix_cell{i} = cropped_CP_Idx_Matrix;
    FC_nN_nm_mat_cell{i}  = cropped_FC_nN_nm;
    
    % For pores cropped near the edge of the image
    % Smaller matrices are cropped if the pores are near the edge of the image.
    % Need to therefore crop a CircleMatrix of the same size, with the 0nm
    % pixel at the same coordinate position as the one designated previously by
    % the user.


    % Find the minimum row/col value and max row/col value. Outside of this
    % range and a smaller matrix will be cropped.
    Crop_coord_max = entireImg_row - (CropSize_pix);
    Crop_coord_min = 1;

    % c1 is col, c2 is row
    if c2>=Crop_coord_min && c2<=Crop_coord_max % if row cropping in right range, centre row is where you'd expect
        crop_row_cent = round(CropSize_pix/2)+1;
    elseif c2<Crop_coord_min
    %    centre_crop_row = c2 +(round(CropSize_pix/2))+1;
        crop_row_cent = c2 +(round(CropSize_pix/2));
    elseif c2>Crop_coord_max
        crop_row_cent = round(CropSize_pix/2)+1;
    end

    % c1 is col, c2 is row
    if c1>=Crop_coord_min && c1<=Crop_coord_max % if row cropping in right range, centre row is where you'd expect
        crop_col_cent = round(CropSize_pix/2)+1;
    elseif c1<Crop_coord_min
    %     centre_crop_col = c1 +(round(CropSize_pix/2))+1;
        crop_col_cent = c1 +(round(CropSize_pix/2));
    elseif c1>Crop_coord_max
        crop_col_cent = round(CropSize_pix/2)+1;
    end

    % The coordinate of the pixel designated by the user to be the central axis
    % of rotation is then required for rotational averaging later. This is
    % stored in 'centres_cropped'.
    centres_cropped_idx(i, :) = [crop_row_cent, crop_col_cent];
    % viscircles uses xy coordinates rather than row, col convention
    centres_cropped_xy(i, :) = [crop_col_cent, crop_row_cent];
    
    % Create a normal CircleMatrix with 0nm at the centre, to be used to
    % crop from for pores near the edge of the image
    Crop_size_matrix = zeros(CropSize_pix+1, CropSize_pix+1);
    [~, CircleMatrix, ~] = NanoMechFuncs.RadialBins(CropSize_nm, Crop_size_matrix, Initial_Radius_nm, Radius_nm, NumberOfBins);
    
    % Pull out the coordinate of the central pixel in CircleMatrix
    [circ_row_cent, circ_col_cent] = find(CircleMatrix == min(CircleMatrix(:)));
    % Get the size of the crops
    [cropped_img_row, cropped_img_col] = size(heightdata_cropped_nm);

    % get number of pixels either side of the central pixel in the cropped image
    bottom_diff_row = crop_row_cent - 1;
    top_diff_row    = cropped_img_row - crop_row_cent;
    bottom_diff_col = crop_col_cent - 1;
    top_diff_col    = cropped_img_col - crop_col_cent;

    % index out a CircleMatrixCrop of the same size as the cropped image, with
    % the central 0nm pixel overlapping with the central axis of rotation, as
    % designated by the user
    CircleMatrixCrop = CircleMatrix(circ_row_cent-bottom_diff_row:circ_row_cent+top_diff_row,...
        circ_col_cent-bottom_diff_col:circ_col_cent+top_diff_col);
    
    % save this circle matrix
    CircleMatrix_cell{i} = CircleMatrixCrop;
    
    % use this potentially shifted circle matrix to get the coordinates for
    % the radial bins, as well as the radial values (which should be the
    % same whether the 0nm is the central pixel or not).
    [Circle_coords, Radial_Value] = NanoMechFuncs.RadialCoordinatesfromCircleMatrices(CircleMatrixCrop, Initial_Radius_nm, Radius_nm, NumberOfBins);

    % Save out the radial values (for plotting later) and Circle_coords (for binning force curves)
    Radial_Values_cell{i} = Radial_Value;
    Circle_coords_cell{i} = Circle_coords;
    
    % Use this coordinate system to radially bin the height, YM, force
    % curves, CP_idx, and count matrices
    
    % First, pre-allocate for radially binned data (rd)
    heightdata_nm_rd_cropped{i} = cell(1, NumberOfBins);
    YM_Matrix_MPa_rd_cropped{i} = cell(1, NumberOfBins);
    FC_nN_nm_mat_rd_cropped{i}  = cell(1, NumberOfBins);
    CP_Idx_Matrix_rd_cropped{i} = cell(1, NumberOfBins);
    FC_count_mat_rd_cropped{i}  = cell(1, NumberOfBins);
    
    for n = 1:length(Circle_coords)
        
        % pre-allocate length of each element in the cell array (will be number of circle coordinates for that radial bin)
        heightdata_nm_rd_cropped{i}{n} = zeros(1, length(Circle_coords(:,end)));
        YM_Matrix_MPa_rd_cropped{i}{n} = zeros(1, length(Circle_coords(:,end)));
        FC_nN_nm_mat_rd_cropped{i}{n}  = cell(1, length(Circle_coords(:,end)));
        CP_Idx_Matrix_rd_cropped{i}{n} = zeros(1, length(Circle_coords(:,end)));
        FC_count_mat_rd_cropped{i}{n}  = zeros(1, length(Circle_coords(:,end)));
        
        for j = 1:length(Circle_coords{n})
            
            % Pull-out the coordinates (idx, not xy)
            coord_row = Circle_coords{n}(j,1);
            coord_col = Circle_coords{n}(j,2);
                        
            % use coordinates to radially bin all the data types
            heightdata_nm_rd_cropped{i}{n}(j) = heightdata_cropped_nm(coord_row, coord_col);
            YM_Matrix_MPa_rd_cropped{i}{n}(j) = cropped_YM(coord_row, coord_col);
            FC_nN_nm_mat_rd_cropped{i}{n}{j}  = cropped_FC_nN_nm{coord_row, coord_col};
            CP_Idx_Matrix_rd_cropped{i}{n}(j) = cropped_CP_Idx_Matrix(coord_row, coord_col);
            FC_count_mat_rd_cropped{i}{n}(j)  = cropped_FC_count_mat(coord_row, coord_col);
            
        end
        
    end
    
    
end


%% Save the cropped height data, contact point indices, and the Young's moduli

% save the original matrices (before cropping)
Height_YM_CP_Cropped_rb.entirematrices.heightdata_nm          = heightdata_nm;
Height_YM_CP_Cropped_rb.entirematrices.YM_Matrix_MPa          = YM_Matrix_MPa;
Height_YM_CP_Cropped_rb.entirematrices.FC_nN_nm_mat           = FC_nN_nm_mat;
Height_YM_CP_Cropped_rb.entirematrices.CP_Idx_Matrix          = CP_Idx_Matrix;
Height_YM_CP_Cropped_rb.entirematrices.FC_count_mat           = FC_count_mat;
Height_YM_CP_Cropped_rb.entirematrices.centres_entire_img_idx = centres_entire_img_idx;
Height_YM_CP_Cropped_rb.entirematrices.ScanSize_nm            = ScanSize_nm;

% save all the cropped matrices
Height_YM_CP_Cropped_rb.matrices.heightdata_nm_cell     = heightdata_nm_cell;
Height_YM_CP_Cropped_rb.matrices.YM_MPa_cell            = YM_MPa_cell;
Height_YM_CP_Cropped_rb.matrices.FC_nN_nm_mat_cell      = FC_nN_nm_mat_cell;
Height_YM_CP_Cropped_rb.matrices.CP_Idx_Matrix_cell     = CP_Idx_Matrix_cell;
Height_YM_CP_Cropped_rb.matrices.FC_count_cell          = FC_count_cell;
Height_YM_CP_Cropped_rb.matrices.CircleMatrix_cell      = CircleMatrix_cell;
Height_YM_CP_Cropped_rb.matrices.centres_cropped_idx    = centres_cropped_idx;
Height_YM_CP_Cropped_rb.matrices.centres_cropped_xy     = centres_cropped_xy;

% save all the radially binned data
Height_YM_CP_Cropped_rb.radiallybinned.Radial_Values_cell     = Radial_Values_cell;
Height_YM_CP_Cropped_rb.radiallybinned.Circle_coords_cell     = Circle_coords_cell;
Height_YM_CP_Cropped_rb.radiallybinned.heightdata_nm_rd_cropped = heightdata_nm_rd_cropped;
Height_YM_CP_Cropped_rb.radiallybinned.YM_Matrix_MPa_rd_cropped = YM_Matrix_MPa_rd_cropped;
Height_YM_CP_Cropped_rb.radiallybinned.FC_nN_nm_mat_rd_cropped  = FC_nN_nm_mat_rd_cropped;
Height_YM_CP_Cropped_rb.radiallybinned.CP_Idx_Matrix_rd_cropped = CP_Idx_Matrix_rd_cropped;
Height_YM_CP_Cropped_rb.radiallybinned.FC_count_mat_rd_cropped  = FC_count_mat_rd_cropped;


FullFileOutput = fullfile(OutputFolder, strcat(GenericFileName, ' - ', ' Height_YM_CP_Cropped_rb', '_', num2str(FileNumb)));
save(strcat(FullFileOutput, '.mat'), 'Height_YM_CP_Cropped_rb');

%% Plot centres and show to user for verification

for i = 1:length(centres_cropped_xy(:,end))
    
    figure();
        
    img = heightdata_nm_cell{i};
    
    [rr, cc] = size(img);
    
    imagesc(img)
    caxis([height_min_nm height_max_nm])
    colorbar;
    g = colorbar;
    ylabel(g,'Height (nm)', 'FontSize', 13)
    set(gca,'ydir','normal');
    set(gca, 'FontSize', 13)
    title('Height data (nm)')
    pbaspect([cc rr 1])
    hold on
    viscircles(centres_cropped_xy(i,:), 1,'EdgeColor','r', 'LineWidth', 0.2);
    
    figure();
        
    ym = YM_MPa_cell{i};
    
    imagesc(ym)
    caxis([0 4])
    colorbar;
    g = colorbar;
    ylabel(g,'YM', 'FontSize', 13)
    set(gca,'ydir','normal');
    set(gca, 'FontSize', 13)
    title('YM')
    pbaspect([cc rr 1])
    hold on
    viscircles(centres_cropped_xy(i,:), 1,'EdgeColor','r', 'LineWidth', 0.2);
    
end
%

radii_1 = ones(length(centres_entire_img_xy(:,end)),1);
% radii_2 = ones(length(centres_entire_img_xy(:,end)),1) .* 2;
radii_3 = ones(length(centres_entire_img_xy(:,end)),1) .* 3;

figure(), subplot(121)
imagesc(heightdata_nm)
set(gca,'ydir','normal');
pbaspect([1 1 1])
hold on
viscircles(centres_entire_img_xy, radii_1,'EdgeColor','r', 'LineWidth', 0.2);
% viscircles(centres_entire_img_xy, radii_2,'EdgeColor','w', 'LineWidth', 0.2);
viscircles(centres_entire_img_xy, radii_3,'EdgeColor','r', 'LineWidth', 0.2);
caxis([height_min_nm height_max_nm])
colorbar
g = colorbar;
ylabel(g,'Height (nm)', 'FontSize', 13)
set(gca,'ydir','normal');
set(gca, 'FontSize', 13)
title('Height (nm)')
hold off

subplot(122), imagesc(YM_Matrix_MPa)
set(gca,'ydir','normal');
pbaspect([1 1 1])
hold on
viscircles(centres_entire_img_xy, radii_1,'EdgeColor','r', 'LineWidth', 0.2);
% viscircles(centres_entire_img_xy, radii_2,'EdgeColor','w', 'LineWidth', 0.2);
viscircles(centres_entire_img_xy, radii_3,'EdgeColor','r', 'LineWidth', 0.2);
caxis([0 4])
colorbar
g = colorbar;
ylabel(g,'YM (MPa)', 'FontSize', 13)
set(gca,'ydir','normal');
set(gca, 'FontSize', 13)
title('YM (MPa)')
hold off

figure(), imagesc(heightdata_nm)
caxis([height_min_nm height_max_nm])
colorbar;
g = colorbar;
ylabel(g,'Height (nm)', 'FontSize', 13)
set(gca,'ydir','normal');
set(gca, 'FontSize', 13)
title('Height data (nm)')
pbaspect([1 1 1])
hold on
viscircles(centres_entire_img_xy, radii_1,'EdgeColor','r', 'LineWidth', 0.2);
% viscircles(centres_entire_img_xy, radii_2,'EdgeColor','w', 'LineWidth', 0.2);
viscircles(centres_entire_img_xy, radii_3,'EdgeColor','r', 'LineWidth', 0.2);



