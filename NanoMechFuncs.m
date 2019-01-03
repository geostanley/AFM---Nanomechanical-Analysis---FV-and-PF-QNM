classdef NanoMechFuncs
    methods(Static)
        
        function [FileNames_spm, FileNames_pfc] = Create_FileNames_Cell_3_spm_pfc(Data_Directory, File_Suffix, File_Nos_spm, File_Nos_pfc, NS_nine_one)

        %%%=== Create_FileNames_Cell_3_spm_pfc ===%%%

        % This function creates a cell containing strings of filenames. Input the
        % directory contatining the desired data files, the file suffix common to
        % all data files, and an array containing the numbers of all file ends.
        % It creates both spm and pfc file names.

        % Also tell the function which version of Nanoscope the spm files are from.
        % If 9.1, input the string 'yes'. If 9.2 or higher, input 'no'.

        % A cell array with all the file names will be returned.

        File_Numbers_spm = cell(1,length(File_Nos_spm));
        for i=1:length(File_Nos_spm);
            if File_Nos_spm(i) < 10
                File_Numbers_spm{i} = strcat('00', num2str(File_Nos_spm(i)));
            elseif File_Nos_spm(i) >=10 && File_Nos_spm(i) < 100
                File_Numbers_spm{i} = strcat('0', num2str(File_Nos_spm(i)));
            else
                File_Numbers_spm{i} = num2str(File_Nos_spm(i));
            end
        end

        File_Numbers_pfc = cell(1,length(File_Nos_pfc));
        for i=1:length(File_Nos_pfc);
            if File_Nos_pfc(i) < 10
                File_Numbers_pfc{i} = strcat('00', num2str(File_Nos_pfc(i)));
            elseif File_Nos_pfc(i) >=10 && File_Nos_pfc(i) < 100
                File_Numbers_pfc{i} = strcat('0', num2str(File_Nos_pfc(i)));
            else
                File_Numbers_pfc{i} = num2str(File_Nos_pfc(i));
            end
        end


        if length(NS_nine_one) == 3

            FileNames_spm = cell(1,length(File_Nos_spm));
            FileNames_pfc = cell(1,length(File_Nos_pfc));

            for i = 1:length(File_Nos_spm)

                FileNames_spm{i} = horzcat(Data_Directory, File_Suffix,'.', File_Numbers_spm{i});
                FileNames_pfc{i} = horzcat(Data_Directory, File_Suffix,'.', File_Numbers_pfc{i}, '.pfc');

            end

        else

            FileNames_spm = cell(1,length(File_Nos_spm));
            FileNames_pfc = cell(1,length(File_Nos_pfc));

            for i = 1:length(File_Nos_spm)

                FileNames_spm{i} = horzcat(Data_Directory, File_Suffix,'.0_00', File_Numbers_spm{i}, '.spm');
                FileNames_pfc{i} = horzcat(Data_Directory, File_Suffix,'.0_00', File_Numbers_pfc{i}, '.pfc');

            end

        end

        end
        
        %% Create_FileNames_Cell_3_any
        
        function [FileNames_any] = Create_FileNames_Cell_3_any(Data_Directory, File_Suffix, File_Nos, FileType)
            
        % This function creates a cell containing strings of filenames. Input the
        % directory contatining the desired data files, the file suffix common to
        % all data files, and an array containing the numbers of all file ends -
        % and the fily type (including the point). It creates a cell array
        % containing all the fullfile names ready to be loaded.

            FileNames_any = cell(1,length(File_Nos));

            for i = 1:length(File_Nos)

                FileNumb = File_Nos(i);

                FileNames_any{i} = fullfile(Data_Directory, strcat(File_Suffix, num2str(FileNumb), FileType));

            end


        end
            
        %% HertzModel
    
        function [YM] = HertzModel(Force_Curve_Corrected, Contact_Index, Indentation_Hertz_Fit)
            
        % This function applies the Hertz model to a force curve. The force curve
        % should be already corrected for the baseline and contact point, and
        % should be in the format nN vs nm (tip-sample separation). The input
        % arguments for the function are: the corrected force curve, entered as an
        % Nx2 array; the index for its contact point along the x axis (tip-sample
        % separation); and the length of indentation, in nm, that is to be
        % modelled. Typically, 20 nm is modelled. The output is Young's Modulus
        % in MPa.

        % Find the point closest to 20 nm before contact point
        Contact_Point_Corrected_x_nm = Force_Curve_Corrected(Contact_Index,1); % CP in nm
        Shifted_Force_Curve_Corrected = abs(Force_Curve_Corrected(:,1) - (Contact_Point_Corrected_x_nm - Indentation_Hertz_Fit));

        % Find the index of the point closest to CP-Identation_Hertz_Fit
        Indent_20_Index_Corrected = find(Shifted_Force_Curve_Corrected == min(Shifted_Force_Curve_Corrected)); 

        % Save this 20nm region as:
        Force_Curve_Hertz_Fit_Corrected = [Force_Curve_Corrected(Indent_20_Index_Corrected:Contact_Index,1),...
            Force_Curve_Corrected(Indent_20_Index_Corrected:Contact_Index,2)];

        % For fitting, save this data as x and y, having inversed the x values so
        % that they range from 0-20nm, rather than -20-0nm. This avoids having sqrt(-x)
        % in the Hertz Model.
        Force_Curve_Hertz_Fit_Corrected_Positive = zeros(length(Force_Curve_Hertz_Fit_Corrected),2);
        Force_Curve_Hertz_Fit_Corrected_Positive(:,1) = Force_Curve_Hertz_Fit_Corrected(:,1) .* -1;
        Force_Curve_Hertz_Fit_Corrected_Positive(:,2) = Force_Curve_Hertz_Fit_Corrected(:,2);
        x = Force_Curve_Hertz_Fit_Corrected_Positive(:,1);
        y = Force_Curve_Hertz_Fit_Corrected_Positive(:,2);

        % Enter Poisson Ratio and Tip Radius.
        PR = 0.5;
        TR = 2;

        % suppress the 'Optimization complete' message
        opts = optimset('Display', 'off');

        % Form of the equation
        Hertz_Model = @(b,x) ((4/3)*(b(1)/(1-PR^2))*sqrt(TR)*(x.^(3/2)));
        % Initial guess at the Young's modulus in GPa (to speed up the least squares regression)
        beta0 = 0.002; % if several parameters use an array
        % Use a non-linear least squares regression function
        YoungsModulus_GPa = nlinfit(x,y,Hertz_Model, beta0, opts);
        % Multiply by 1000 to give it in MPa
        YoungsModulus_MPa = YoungsModulus_GPa*1000;
        % Ignore any imaginary components (produced by taking the square root of a negative number)
        YoungsModulus_MPa = real(YoungsModulus_MPa);

        YM = YoungsModulus_MPa;

        %%%===========================End Hertz=================================%%%

        end
        
        %% InputCentres_Colourmap
        
        function Centres = InputCentres_Colourmap(Image, ColourMap)
            
        % This function is used to assign the coordinates for the central axes of
        % rotation for the NPCs. The input argument is the AFM image in greyscale.
        % A figure is then produced showing the image, and the user must manually
        % click on the picture to assign the central axes of rotation for each NPC.
        % The output is an Nx2 array of the coordinates, to the nearest pixel, of
        % the selected centres.


        % for saturating colour scale, first take histogram of height data, binning
        % into 1nm bins

        heightdata=Image;
        heightdata_array = heightdata(:);
        histedges = [floor(min(heightdata_array)):ceil(max(heightdata_array))];
        hist_counts = histcounts(heightdata_array, histedges);
        hist_counts_sum = sum(hist_counts(:));

        % create an array with aggregate counts based on bins
        aggregate_hist_counts = zeros(size(hist_counts));
        aggregate_hist_counts(1) = hist_counts(1);
        for i=1:length(hist_counts)-1
            aggregate_hist_counts(i+1) = aggregate_hist_counts(i) + hist_counts(i+1);
        end

        % find bin which contains 98% of aggregate counts (or closest to)
        hist_counts_sum_ninety_five = hist_counts_sum*0.98;
        aggregate_hist_counts_abs = abs(aggregate_hist_counts - round(hist_counts_sum_ninety_five));
        cut_off_idx = aggregate_hist_counts_abs == min(aggregate_hist_counts_abs);

        % find bin edges (in nm) for this 98% of aggregate count, and define the
        % colour limits in nm
        colour_lim_floor = histedges(1);
        colour_lim_ceil = histedges(cut_off_idx);
        clims = [colour_lim_floor colour_lim_ceil];

        
        figure();
        imagesc(Image, clims);
        pbaspect([1 1 1])
        colormap(ColourMap);
        caxis(clims);
        set(gca,'ydir','normal');
        xlabel('Select centres and press enter when finished')
        set (gcf, 'WindowButtonMotionFcn', @mouseMove);
        [c, r] = ginput;
        r = round(r);
        c = round(c);
        Centres = [c r];


        end
        
        %% RadialBins
        
        function [Circle_coords, Circle_Matrix, Radial_Value] = RadialBins(Crop_Size, Matrix_Of_Same_Size, Initial_Radius, Radius, NumberOfBins)
            
        %%%========================Create Radial Bins===========================%%%

        % This function gives a cell with the coordinates of element positions for
        % the different radial bins. The input arguments are: the size of the image
        % in nm; the size of the matrix in pixels (this is input as any matrix of
        % the correct size); the radius for the first bin, in nm; the radius for
        % all subsequent bins; and finally, the number of bins in total. The
        % outputs are: the coordinates of the bins as a cell; the matrix with the
        % nm value corresponding to each pixel; and an array containing the average
        % radial value, in nm, for each bin.

        % Rows and columns of a cropped image
        [rr, cc] = size(Matrix_Of_Same_Size);

        Ones_Matrix = ones(rr,cc);

        % nm/pixels for 130nm crop
        Ratio_nm_per_pixs = Crop_Size/rr;

        % Pre-allocate cell array for the coordinates of the pixels of the circles.
        % The number of cells will be the number of bins.
        Circle_coords      = cell(1,NumberOfBins);
        Mask_Circle        = cell(1,NumberOfBins);
        Mask_Store         = cell(1,NumberOfBins);  
        Mask_Circle_Pixels = cell(1,NumberOfBins);
        Radial_Value       = zeros(1,NumberOfBins);

        % Manually create the first bin (0-6nm radius)
        min_boundary = 0;
        max_boundary = Initial_Radius;
        [Mesh_x, Mesh_y] = meshgrid(-((rr/2)-0.5):((cc/2)-0.5));
        % C is a matrix whose central element = 0, and all elements leading away
        % increase linearly in a circular pattern.
        C = sqrt(((Mesh_x).^2) + ((Mesh_y).^2));
        % Circle matrix converts these values using the scalar Ratio_nm_per_pixels
        % to convert the matrix to reflect the real nm values.
        Circle_Matrix = C*Ratio_nm_per_pixs;
        % Index the nm values desired.
        Mask_Circle{1} = Circle_Matrix>=min_boundary & Circle_Matrix<max_boundary;
        Mask_Circle_Pixels{1} = Circle_Matrix(Circle_Matrix >=...
                min_boundary & Circle_Matrix < max_boundary);
        % Use the mask to delete all other pixel values.
        Masked_Pore = Ones_Matrix.*Mask_Circle{1};

        [Coords_rr, Coords_cc] = find(Masked_Pore>0);
        Circle_coords{1} = [Coords_rr, Coords_cc];


        Mask_Store{1} = zeros(rr,cc);
        for j=1:length(Coords_rr)
            Mask_Store{1}(Coords_rr(j),Coords_cc(j))=1;
        end

        for i=1:NumberOfBins-1
            Min_Boundary = Initial_Radius + ((i-1)*Radius);
            Max_Boundary = Initial_Radius + (i*Radius);
            Mask_Circle{i+1} = Circle_Matrix >=...
                Min_Boundary & Circle_Matrix < Max_Boundary;
            Mask_Circle_Pixels{i+1} = Circle_Matrix(Circle_Matrix >=...
                Min_Boundary & Circle_Matrix < Max_Boundary);
            Masked_Ones = Ones_Matrix.*Mask_Circle{i+1};
            [Coords_rr, Coords_cc] = find(Masked_Ones>0);
            Circle_coords{i+1} = [Coords_rr, Coords_cc];
            Mask_Store{i+1} = zeros(rr,cc);
            for j=1:length(Circle_coords{i+1})
                Mask_Store{i+1}(Coords_rr(j),Coords_cc(j))=1;
            end
        end


        % Need to save values of pixels in each bin, radial values of circle
        % matrix, so can find average radial value for the radii array.
        for i=1:NumberOfBins
             Radial_Value(i) = mean(Mask_Circle_Pixels{i});
        end

        % overwrite first value and set to 0nm.
        Radial_Value(1) = 0;

        %%%=====================End Create Radial Bins===========================%%%

        end
        
        
        %% RadialCoordinatesfromCircleMatrices
        
        function [Circle_coords, Radial_Value] = RadialCoordinatesfromCircleMatrices(Circle_Matrix, Initial_Radius, Radius, NumberOfBins)
            
        %%%=== RadialCoordinatesfromCircleMatrices ===%%%

        % This function is a variation of the RadialBins function. It takes in a
        % CircleMatrix of any size with a 0nm pixel not necessarily at the centre,
        % and returns the coordinates for each radial bin. The limits of the radial
        % bins are specified by the user.

        Mask_Circle        = cell(1,NumberOfBins);
        Mask_Store         = cell(1,NumberOfBins);  
        Mask_Circle_Pixels = cell(1,NumberOfBins);
        Radial_Value       = zeros(1,NumberOfBins);
        Circle_coords      = cell(1,NumberOfBins);

        [rr, cc] = size(Circle_Matrix);
        Ones_Matrix = ones(rr,cc);

        % Manually create the first bin (0-6nm radius)
        min_boundary = 0;
        max_boundary = Initial_Radius;

        % Index the nm values desired.
        Mask_Circle{1} = Circle_Matrix>=min_boundary & Circle_Matrix<max_boundary;
        Mask_Circle_Pixels{1} = Circle_Matrix(Circle_Matrix >=...
                min_boundary & Circle_Matrix < max_boundary);
        % Use the mask to delete all other pixel values.
        Masked_Pore = Ones_Matrix.*Mask_Circle{1};

        [Coords_rr, Coords_cc] = find(Masked_Pore>0);
        Circle_coords{1} = [Coords_rr, Coords_cc];

        Mask_Store{1} = zeros(rr,cc);
        for j=1:length(Coords_rr)
            Mask_Store{1}(Coords_rr(j),Coords_cc(j))=1;
        end

        for i=1:NumberOfBins-1
            Min_Boundary = Initial_Radius + ((i-1)*Radius);
            Max_Boundary = Initial_Radius + (i*Radius);
            Mask_Circle{i+1} = Circle_Matrix >=...
                Min_Boundary & Circle_Matrix < Max_Boundary;
            Mask_Circle_Pixels{i+1} = Circle_Matrix(Circle_Matrix >=...
                Min_Boundary & Circle_Matrix < Max_Boundary);
            Masked_Ones = Ones_Matrix.*Mask_Circle{i+1};
            [Coords_rr, Coords_cc] = find(Masked_Ones>0);
            Circle_coords{i+1} = [Coords_rr, Coords_cc];
            Mask_Store{i+1} = zeros(rr,cc);
            for j=1:length(Circle_coords{i+1})
                Mask_Store{i+1}(Coords_rr(j),Coords_cc(j))=1;
            end
        end

        % Need to save values of pixels in each bin, radial values of circle
        % matrix, so can find average radial value for the radii array.
        for i=1:NumberOfBins
             Radial_Value(i) = mean(Mask_Circle_Pixels{i});
        end

        % overwrite the first bin and make 0nm
        Radial_Value(1) = 0;    

        end
        
        %% Average_forcecurve 

        function[count_beforeThresh, Averaged_plot_beforeThresh, count_final, Averaged_plot] = Average_forcecurves(fcs, binwidth, threshold)
            
        % This function takes in a cell array of force curves and averages them. It
        % is assumed that all force curves have already been aligned (i.e. by the
        % contact point). After averaging the force curves, it further indexes out
        % only the data points that have contributed to the averaging process a
        % significant amount. It does this by keeping count of how many times a
        % data point is placed in a bin (this counting process is anyway required
        % for the averaging of the force curve), and then dividing by the number of
        % force curves. Then, only data points greater than the threshold value are
        % accepted. E.g., if the threshold value is 0.7, only data points for which
        % at least 70% of the force curves have contributed, will be accepted. This
        % removes non-linear behaviour at both ends of the averaged force curve.

        % The user must also input a sensible binwidth. For averaging force cueves
        % with a ramp size of 100nm, I use a binwidth of 1 (nm).

        % The outputs are the counts and average plots of the force curves before
        % indexing using the threshold value, and after.

        numbfcs = length(fcs);

        % go through all the force curves and search for the minimum and maximum x
        % values. This is to create the bin edges for binning the data later.
        xplotmin_array = zeros(size(fcs));
        xplotmax_array = zeros(size(fcs));

        for i=1:length(fcs)

            xdata = fcs{i}(:,1);

            xplotmin = min(xdata);
            xplotmax = max(xdata);

            xplotmin_array(i) = xplotmin;
            xplotmax_array(i) = xplotmax;

        end

        min_binedge = floor(min(xplotmin_array));
        max_binedge = ceil(max(xplotmax_array));

        % bin the data and keep a count of how many data points placed into each
        % bin
        xbins = min_binedge:binwidth:max_binedge;
        accumulative_sum_y = zeros(1, length(xbins));
        count              = zeros(1, length(xbins));
        for j=1:length(fcs)

            xdata_temp = fcs{j}(:,1);
            ydata_temp = fcs{j}(:,2);

                % put the data points into the correct bins, and keep a count of them. 
                for i=1:length(xdata_temp)
                    closest_x                   = abs(xbins - xdata_temp(i));
                    pos_idx                     = find(closest_x == min(closest_x));
                    accumulative_sum_y(pos_idx) = accumulative_sum_y(pos_idx) + ydata_temp(i);
                    count(pos_idx)              = count(pos_idx) + 1;
                end    

        end

        % any points divided by 0 will appear as nan. To avoid this hapenning, we
        % search in the count index for points greater than 1, and only keep these
        idx_keep = find(count>0);
        accumulative_sum_y_beforeThresh = accumulative_sum_y(idx_keep);
        xbins_beforeThresh              = xbins(idx_keep);
        count_beforeThresh              = count(idx_keep);

        % average the data points and save out the new array
        averaged_ydata_beforeThresh = accumulative_sum_y_beforeThresh./count_beforeThresh;
        Averaged_plot_beforeThresh = [xbins_beforeThresh(:), averaged_ydata_beforeThresh(:)];

        % Now count how many data points have contributed to each bin. If less than a threshold, remove.
        count_keep_norm = count_beforeThresh ./ numbfcs;
        count_final_idx = find(count_keep_norm > threshold);

        accum_y_final = accumulative_sum_y_beforeThresh(count_final_idx);
        xbins_final   = xbins_beforeThresh(count_final_idx);
        count_final   = count_beforeThresh(count_final_idx);

        % average the data points and save out the new array
        averaged_ydata_final = accum_y_final./count_final;
        Averaged_plot = [xbins_final(:), averaged_ydata_final(:)];

        end
        
        %% Average_forcecurves_maximum_indentation
    
        function[FC_ave_max_indentation] = Average_forcecurves_maximum_indentation(fcs, binwidth)
            
        % This function takes in a cell array of force curves and averages them.
        % This is done by simply aligning the force curves by their maximum
        % indentation. It is assumed that their baselines have already been
        % aligned.

        % The user must input a sensible binwidth. For averaging force curves
        % with a ramp size of 100nm, I use a binwidth of 1 (nm).

        % The output is the averaged plots of the force curves.

        numbfcs = length(fcs);

        % Align all force curves on zero nm, and store the length of each one

        fcs_xadjusted     = cell(1, numbfcs);
        xplotmax_adjusted = zeros(1, numbfcs);

        for i = 1:numbfcs

            xplot = fcs{i}(:,1);
            yplot = fcs{i}(:,2);

            xplotmin = min(xplot);

            xplot_adjusted     = xplot - xplotmin;
            xplot_adjusted_max = max(xplot_adjusted);

            xplotmax_adjusted(i) = xplot_adjusted_max;

            fcs_xadjusted{i} = [xplot_adjusted, yplot];

        end

        max_binedge = round(max(xplotmax_adjusted));

        % bin the data and keep a count of how many data points placed into each
        % bin
        xbins = 0:binwidth:max_binedge;

        accumulative_sum_y = zeros(1, length(xbins));
        count              = zeros(1, length(xbins));

        for j=1:length(fcs)

            xdata_temp = fcs_xadjusted{j}(:,1);
            ydata_temp = fcs_xadjusted{j}(:,2);

                % put the data points into the correct bins, and keep a count of them. 
                for i=1:length(xdata_temp)
                    closest_x                   = abs(xbins - xdata_temp(i));
                    pos_idx                     = find(closest_x == min(closest_x));
                    accumulative_sum_y(pos_idx) = accumulative_sum_y(pos_idx) + ydata_temp(i);
                    count(pos_idx)              = count(pos_idx) + 1;
                end    

        end

        % average the data points and save out the new array
        averaged_ydata = accumulative_sum_y./count;
        FC_ave_max_indentation = [xbins(:), averaged_ydata(:)];

        end
        
        %% DerivativeForceCurve_2

        function[Matlab_Stiffness_Curve, Stiffness_Curve] = DerivativeForceCurve_2(ForceCurve, Invert_Matlab_Derivative)
            
        % Enter a cell array of force curves, and a cell array of their derivatives
        % is returned.

        % The results from my own derivative code, and from Matlab's own function,
        % give exactly the same results.

        % Stiffness curves from the averaged force curves, i.e., take the
        % derivative of each force curve.

        % Preallocate the Gradients array
        Stiffness_Curve        = cell(1,length(ForceCurve));
        Matlab_Stiffness_Curve = cell(1,length(ForceCurve));

        for j=1:length(Stiffness_Curve)

            if isempty(ForceCurve{j}) == 0

                Stiffness_Curve{j}            = zeros((length(ForceCurve{j})-1),2);
                Matlab_Stiffness_Curve_x_diff = zeros((length(ForceCurve{j})-1),1);
            %     Matlab_Stiffness_Curve_y      = zeros((length(ForceCurve{j})),1);
            %     Matlab_Stiffness_Curve_y_diff = zeros((length(ForceCurve{j})-1),1);

                for i=1:length(Stiffness_Curve{j})

                    Matlab_Stiffness_Curve_x_diff(i) = (ForceCurve{j}(i,1) + ForceCurve{j}(i+1,1))/2;
                    Stiffness_Curve{j}(i,1)          = (ForceCurve{j}(i,1) + ForceCurve{j}(i+1,1))/2;
                    Stiffness_Curve{j}(i,2)          = (ForceCurve{j}(i,2) - ForceCurve{j}(i+1,2));
                end


                Matlab_Stiffness_Curve_y      = ForceCurve{j}(:,2);
                if length(Invert_Matlab_Derivative) == 3
                    Matlab_Stiffness_Curve_y = Matlab_Stiffness_Curve_y * (-1);
                else
                    continue
                end

                Matlab_Stiffness_Curve_y_diff = diff(Matlab_Stiffness_Curve_y);
                Matlab_Stiffness_Curve{j} = [Matlab_Stiffness_Curve_x_diff, Matlab_Stiffness_Curve_y_diff];

            end

        end

        end
         
    
        
    end  
end