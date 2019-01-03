%%%=== NM_Rotate_Average ===%%%
 
% This is script 4/5 for the nanomechanical analysis procedure.

% This script loads all the concatonated data from the data structure:
% GenericFileName - Height_YM_CP_Cropped_rb_concatonated.mat

% It then averages all the data based on their radial distance from the
% centre of the pore - ready for plotting in the next (and final) script.

% It also saves out the indentation values from every force curve for 
% plotting of the true height. This should have been done in the first script,
% but I forgot, and I've already processed too much data, so I'm throwing it 
% in here like the cowboy I always knew I could be.

%% Enter load and save directories and file names etc

clear variables
close all
clc

display('NM_Rotate_Average')

% can decide whether to save out results in data structure. This is because
% the results can be loaded into separate figure plotting scripts which
% render journal quality plots
save_out_data_structure = 0;

%%%%%%%=== File to be loaded === %%%%%%%%%%%%%%%%%
GenericFileName = '2kHz_test_cyto';

%%%%%%%=== Data structure directory
LoadFolder = ['Z:\Users\George\Documents\PhD\Data\'...
    'Nanomechanical_Outputs_Hertz_YM_CP\Test'];

%%%%%%%=== Output folder
OutputFolder = ['Z:\Users\George\Documents\PhD\Data\'...
    'Nanomechanical_Outputs_Hertz_YM_CP\Test'];
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%=== DO NOT CHANGE ===%%%%%%%%%%%%
DataStructureName = ' - Height_YM_CP_Cropped_rb_concatonated';
LoadFileName = strcat(GenericFileName, DataStructureName, '.mat');

FullFileName = fullfile(LoadFolder, LoadFileName);

%%%=== Parameters for averaging force curves ===%%%
binwidth_nm = 1;   % decide on a sensible bin width for placing the force curve data into bins
threshold   = 0.8; % only keep the indentation values for which x% of fcs contributed (usually I go for 80%, i.e. 0.8)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HertzIndentation_nm = 20;

% Make nice colours for plotting
N = 4;
C = linspecer(N);
c1 = C(1,:);
c2 = C(2,:);
c3 = C(3,:);
c4 = C(4,:);

%% Load data structure

display('Loading data structure...')
load(FullFileName)

% pull-out the radially binned heightdata_nm (hd), Young's moduli (ym),
% force curves (nN/nm, corrected; fc), and count.
hd    = radiallybinned.heightdata_nm_rd_cropped_cat;
ym    = radiallybinned.YM_Matrix_MPa_rd_cropped_cat;
fc    = radiallybinned.FC_nN_nm_mat_rd_cropped_cat;
cc    = radiallybinned.FC_count_mat_rd_cropped_cat;
% get the radial values.
rv    = radiallybinned.Radial_Values_cell_cat{1};

numbpores = length(hd);

%% First, find indentation values
% this is because will use indentation values to bin force curves for which
% contact point found too near either extreme of force curve (which would
% suggest a very strange forve curve).

display('Finding indentation values and contact point indices...')

iv = cell(size(fc));

for i = 1:length(fc)
        
   for j = 1:length(fc{i})
              
       for n = 1:length(fc{i}{j})
           
           fc_current = fc{i}{j}{n};
           fc_x       = fc_current(:,1);
           ind_val    = abs(min(fc_x));
           
           iv{i}{j}(n) = ind_val;

       end
   end
    
end

% go through the indentation values, where outside the range, bin the force
% curve and change the counter

cc_iv = cc;
fc_iv = fc;
ym_iv = ym;

for i =1:length(iv)
    for j = 1:length(iv{i})
        for n = 1:length(iv{i}{j})
            
            ind_value = iv{i}{j}(n);
            
            if ind_value<=5 || ind_value>=75 || isnan(ind_value) == 1 % if indentation tiny or huge, or if the contact point was never found
                
                fc_iv{i}{j}{n} = nan; % if true, bin force curve
                ym_iv{i}{j}(n) = 0; % bin concomitant ym value
                cc_iv{i}{j}(n) = 0; % and, keep track of binning of force curves in count
                
            end
        end
    end
end
                
                
% Now, go through the YM data and index for 0.1MPa<=YM<10mPa
% only want to keep the YM values and their concomitant force curves for
% values which lie within a reasonable range. This should hopefully bin
% force curves that have interacted directly with the glass (thus giving a
% very high YM), and force curves which never really indented into the
% material.

% Therefore, go through YM values, and where out of range, change the YM
% value to 0, create a new count cell to keep track, and remove the
% concomitant force curves (change them to nans).

cc_iv_ymidx = cc_iv;
fc_iv_ymidx = fc_iv;
ym_iv_ymidx = ym_iv;

for i = 1:length(ym_iv)
    for j = 1:length(ym_iv{i})
        
        for n = 1:length(ym_iv{i}{j})
            
            ym_value  = ym_iv{i}{j}(n);
            cc_value  = cc_iv{i}{j}(n);
            ind_value = iv{i}{j}(n);
            
            if ym_value<=0.1 || ym_value>=10 
                ym_value = 0;
                cc_value = 0;               
                fc_iv_ymidx{i}{j}{n} = nan;
            end
            
            cc_iv_ymidx{i}{j}(n) = cc_value;
            ym_iv_ymidx{i}{j}(n) = ym_value;
            
        end
    end
end


%% Re-order force curves into radial bins (this will be done later for height data and YM values)
% this is slightly more complex than for height or YM, as some force curves
% have been chucked. This section goes through all the data, only storing
% the force curves that exist, leaving no gaps.

% the same is done for the indentation values

fc_rb    = cell(length(rv), 1);
iv_rb    = cell(length(rv), 1);

length_memory_fc = zeros(length(rv), 1);

for i = 1:numbpores % number of pores
    for j = 1:length(fc_iv_ymidx{i}) % number of bins       
                
        count = cc_iv_ymidx{i}{j};
        numb_fcs = sum(count(:)); % this is number of force curves in this bin of this pore
                
        % as cropped matrices can vary in size, so can the length of the
        % same radial bins between pores. Therefore must keep track of how
        % many pixels are in each radial bin to ensure we don't add in
        % extra zeros which will affect the averaging process.
        previous_length = length_memory_fc(j);
        length_memory_fc(j) = previous_length + numb_fcs;
        
        current_length = length_memory_fc(j);
        
        % as some force curves have been removed and replaced with an nan,
        % we need an fc_nan_counter to track them, and adjust the indexing
        % so that we only pull-out the real force curves
        fc_nan_counter = 0;
        
        if numb_fcs > 0 % if there are fcs in the bin
            
            for n = 1:length(count); % for length of pixels in bin (binned or not)
                
                if count(n) == 1 % and if the force curve exists           

                    fc_current = fc_iv_ymidx{i}{j}(n);
                    fc_rb{j}(1, (current_length - numb_fcs) + (n - fc_nan_counter)) = fc_current; % hd re-organised into radial bin positions
                    
                    iv_value = iv{i}{j}(n);
                    iv_rb{j}(1, (current_length - numb_fcs) + (n - fc_nan_counter)) = iv_value;

                
                else % if no force curve exists, the counter increases by 1, so that in the next iteration fc_current will be placed in idx n-1, thereby not skipping a position                 
                    fc_nan_counter = fc_nan_counter + 1;
                    continue
                end
                
            end
            
        end

    end
           
end

%% Re-order each data set into Nx1 cell arrays, with each row a radial bin

hd_rb = cell(length(rv), 1);
ym_rb = cell(length(rv), 1);
cc_rb = cell(length(rv), 1);

length_memory = zeros(length(rv), 1);

for i = 1:numbpores % number of pores
    for j = 1:length(hd{i}) % number of bins

        numb_pixs = length(hd{i}{j}); % number of pixels in bin (this will be same for ym, cc, and fc)
        
        % as cropped matrices can vary in size, so can the length of the
        % same radial bins between pores. Therefore must keep track of how
        % many pixels are in each radial bin to ensure we don't add in
        % extra zeros which will affect the averaging process.
        previous_length = length_memory(j);
        length_memory(j) = previous_length + numb_pixs;
        
        current_length = length_memory(j);
        
        for n = 1:numb_pixs; % for length of pixels in bin
            
            hd_value = hd{i}{j}(n);
            ym_value = ym_iv_ymidx{i}{j}(n);
            cc_value = cc_iv_ymidx{i}{j}(n);
            hd_rb{j}(1, (current_length - numb_pixs) + n) = hd_value; % hd re-organised into radial bin positions
            ym_rb{j}(1, (current_length - numb_pixs) + n) = ym_value;
            cc_rb{j}(1, (current_length - numb_pixs) + n) = cc_value;
            
        end
    end
end



%% Average the height data
% all height pixels available, so averaging is simple

hd_rb_ave = zeros(1, length(hd_rb));
for i = 1:length(hd_rb)
    
    hd_bin_all = hd_rb{i}; 
    hd_bin_ave = mean(hd_bin_all(:));
    hd_rb_ave(i) = hd_bin_ave;
        
end

%% Average YM values (must be done separately to hd as some force curves may have been removed)
% for averaging the ym values, must use the count array as some force 
% curves were removed

ym_rb_ave         = zeros(1, length(ym_rb));
cc_rb_length      = zeros(1, length(ym_rb));
cc_rb_sum         = zeros(1, length(ym_rb));
percent_rb_bin    = zeros(1, length(ym_rb));

cp_rb_ave         = zeros(1, length(ym_rb));

for i = 1:length(cc_rb)
    
    cc_len = length(cc_rb{i});
    cc_sum = sum(cc_rb{i}(:));
    
    cc_rb_length(i) = cc_len;
    cc_rb_sum(i)    = cc_sum;
    percent_rb_bin(i) = 100 - (cc_sum/cc_len)*100;
    
    ym_bin_all = ym_rb{i};
    ym_bin_sum = sum(ym_bin_all(:));
    ym_rb_ave(i)  = ym_bin_sum/cc_sum;
    
    
end

%% Calculate the standard deviation from the radially binned ym and hd values

% first, create a new cell of ym values, removing the zeros
ym_rb_std_bin = cell(1, length(ym_rb));
length_memory = zeros(1, length(ym_rb));

for i = 1:length(ym_rb)
    
    ym_nan_counter = 0;
    
    for j = 1:length(ym_rb{i})
           
        cc_score = cc_rb{i}(j);
        ym_value = ym_rb{i}(j);
        
        if cc_score == 1;
            
            ym_rb_std_bin{i}(1, (j - ym_nan_counter)) = ym_value;
            
        else
            
            ym_nan_counter = ym_nan_counter + 1;
            
        end
    end
end

% now, calculate std of hd and ym values
ym_rb_std = zeros(1, length(ym_rb_std_bin));
hd_rb_std = zeros(1, length(ym_rb_std_bin));

for i = 1:length(ym_rb_std)
    
    ym_rb_array  = ym_rb_std_bin{i}(:);    
    ym_std       = std(ym_rb_array);    
    ym_rb_std(i) = ym_std;
    
    hd_rb_array  = hd_rb{i}(:);    
    hd_std       = std(hd_rb_array);    
    hd_rb_std(i) = hd_std;
    
end

%% mirror and concatonate for plotting

ym_rb_std_mirror   = fliplr(ym_rb_std(2:end));
ym_rb_std_plot     = [ym_rb_std_mirror, ym_rb_std];
ym_rb_halfstd_plot = ym_rb_std_plot .* 0.5;

hd_rb_std_mirror   = fliplr(hd_rb_std(2:end));
hd_rb_std_plot     = [hd_rb_std_mirror, hd_rb_std];
hd_rb_halfstd_plot = hd_rb_std_plot .* 0.5;

rv_mirror = (-1)*fliplr(rv(2:end));
hd_mirror = fliplr(hd_rb_ave(2:end));
ym_mirror = fliplr(ym_rb_ave(2:end));

rv_plot = [rv_mirror, rv];
hd_plot = [hd_mirror, hd_rb_ave];
ym_plot = [ym_mirror, ym_rb_ave];

% plot average height with error bars (1st standard deviation)
figure(), errorbar(rv_plot, hd_plot, hd_rb_halfstd_plot, 'Color', c2, 'LineWidth', 1)
xlabel('Radial position (nm)', 'FontSize', 13)
ylabel('Height (nm)', 'FontSize', 13)
title(['Height (averaged)', ' (n=', num2str(numbpores), ')'], 'FontSize', 14)
xlim([-100 100])
ylim([0 60])
pbaspect([2 1 1])

% plot average E_eff with error bars (1st standard deviation)
figure(), errorbar(rv_plot, ym_plot, ym_rb_halfstd_plot, 'Color', c1, 'LineWidth', 1)
xlabel('Radial position (nm)', 'FontSize', 13)
ylabel('E_{eff} (MPa)', 'FontSize', 13)
title(['E_{eff} (averaged)', ' (n=', num2str(numbpores), ')'], 'FontSize', 14)
xlim([-100 100])
ylim([0 6])
pbaspect([2 1 1])

%% go through indentation values

iv_rb_min = zeros(1, length(iv_rb));
iv_rb_max = zeros(1, length(iv_rb));
iv_rb_ave = zeros(1, length(iv_rb)); % this is for getting the averaged contact point of force curves aligned by maximum indentation

for i = 1:length(iv_rb)
    
    iv_array = iv_rb{i}(:);
    
    iv_array_length = length(iv_array);
    iv_array_sum    = sum(iv_array(:));
    iv_ave = iv_array_sum/iv_array_length;
    
    iv_rb_ave(i) = iv_ave;
   
    iv_rb_min(i) = min(iv_array);
    iv_rb_max(i) = max(iv_array);
    
end

iv_rb_min_overall = floor(min(iv_rb_min(:)));
iv_rb_max_overall = ceil(max(iv_rb_max(:)));

cutoff_iv_nm_rb = zeros(1, length(iv_rb));

histedges = iv_rb_min_overall:iv_rb_max_overall;

for i = 1:length(iv_rb)
    
    iv_array = iv_rb{i}(:);

    hist_counts = histcounts(iv_array, histedges);
    hist_counts_sum = sum(hist_counts(:));

    % create an array with aggregate counts based on bins
    aggregate_hist_counts = zeros(size(hist_counts));
    aggregate_hist_counts(1) = hist_counts(1);
    for j=1:length(hist_counts)-1
        aggregate_hist_counts(j+1) = aggregate_hist_counts(j) + hist_counts(j+1);
    end

    % find bin which contains 98% of aggregate counts (or closest to)
    hist_counts_sum_thresh = hist_counts_sum * (1-threshold);
    aggregate_hist_counts_abs = abs(aggregate_hist_counts - hist_counts_sum_thresh);
    cut_off_idx = find(aggregate_hist_counts_abs == min(aggregate_hist_counts_abs), 1);

    cutoff_iv_nm = histedges(cut_off_idx);

    cutoff_iv_nm_rb(i) = cutoff_iv_nm;
    
end

%% Average force curves, save indentation values, calculate true height

display('Averaging force curves, getting indentation values, and calculating the true height...')

% pre-allocate
fcs_nN_ave    = cell(1, length(rv));
fcs_pN_ave    = cell(1, length(rv));
fcs_nN_ave_mi = cell(1, length(rv));
fcs_pN_ave_mi = cell(1, length(rv));
ind_rb_ave    = zeros(size(fcs_nN_ave));
th_rb_ave     = zeros(size(fcs_nN_ave));

for i = 1:length(fcs_nN_ave)
    
    % average fcs by contact point
    fcs_rb_to_ave = fc_rb{i};
    [~, fc_averaged, ~, ~] = NanoMechFuncs.Average_forcecurves(fcs_rb_to_ave, binwidth_nm, threshold);
    
    cutoff_iv_nm  = (-1) * cutoff_iv_nm_rb(i);
    fc_averaged_x = fc_averaged(:,1);
    
    idx          = find(fc_averaged_x == cutoff_iv_nm);  
    idx_baseline = find(fc_averaged_x == (fc_averaged_x(end)-20));
    
    if isempty(idx) == 1
        idx = 1;
    end
    
    fc_averaged_cropped = fc_averaged(idx:idx_baseline, :);
    
    fcs_nN_ave{i} = fc_averaged_cropped;
    % also convert from nN into pN
    fc_x = fc_averaged_cropped(:,1);
    fc_y = fc_averaged_cropped(:,2);
    fc_y_pN = fc_y .* 1000;
    fcs_pN_ave{i} = [fc_x, fc_y_pN];
    
    % averaged fcs by maximum indentation
    [fc_averaged_mi] = NanoMechFuncs.Average_forcecurves_maximum_indentation(fcs_rb_to_ave, binwidth_nm);
    fcs_nN_ave_mi{i} = fc_averaged_mi;
    fc_x_mi = fc_averaged_mi(:,1);
    fc_y_mi = fc_averaged_mi(:,2);
    fc_y_pN_mi = fc_y_mi .* 1000;
    fcs_pN_ave_mi{i} = [fc_x_mi, fc_y_pN_mi];


    indentation   = cutoff_iv_nm_rb(i);
    ind_rb_ave(i) = indentation;
    
    trueheight   = hd_rb_ave(i) + indentation;
    th_rb_ave(i) = trueheight;
    
    
end

%% Find averaged contact point index for force curves aligned by max indentation, and apply Hertz model

cp_idx_ave_mi = zeros(size(iv_rb_ave));

for i = 1:length(fcs_nN_ave_mi)
    
    fc_nN_mi = fcs_nN_ave_mi{i};
    fc_nN_mi_x = fc_nN_mi(:,1);
    
    ind_val = round(iv_rb_ave(i));
    
    fc_nN_mi_x_abs = abs(fc_nN_mi_x - ind_val);
    cp_idx_val = find(fc_nN_mi_x_abs == min(fc_nN_mi_x_abs), 1);
    
    cp_idx_ave_mi(i) = cp_idx_val;
    
end

%% Show force curves (aligned by max indentation) and found contact points

figure(), hold on
for i = 1:length(fcs_pN_ave_mi)
    
    fc_x_mi = fcs_pN_ave_mi{i}(:,1);
    fc_y_mi = fcs_pN_ave_mi{i}(:,2);
    
    cp_idx = cp_idx_ave_mi(i);
    
    plot(fc_x_mi, fc_y_mi, 'k')
    scatter(fc_x_mi(cp_idx), fc_y_mi(cp_idx), 'g', 'x')
end
plot(fcs_pN_ave_mi{1}(:,1), fcs_pN_ave_mi{1}(:,2), 'r', 'LineWidth', 2)
scatter(fcs_pN_ave_mi{1}(cp_idx_ave_mi(1), 1), fcs_pN_ave_mi{1}(cp_idx_ave_mi(1), 2), 'r', 'x', 'LineWidth', 2)
xlabel('Z (nm)', 'FontSize', 13)
ylabel('Force (pN)', 'FontSize', 13)
title('Force curves (MI) - found contact points', 'FontSize', 14)
pbaspect([2 1 1])
hold off

%% Correct fcs and remove same amount of indentation values as was done for fcs averaged by cp

% correct fcs
fcs_nN_ave_mi_corrected = cell(size(fcs_nN_ave_mi));

for i = 1:length(fcs_nN_ave_mi)
    
    fc_nN_mi = fcs_nN_ave_mi{i};
        
    fc_nN_mi_x = fc_nN_mi(:,1);
    fc_nN_mi_y = fc_nN_mi(:,2);
    
    cpidx = cp_idx_ave_mi(i);
    cp_nm = fc_nN_mi_x(cpidx);
    
    fc_nN_mi_x_corrected = fc_nN_mi_x - cp_nm;
    
    fc_nN_mi_corrected = [fc_nN_mi_x_corrected, fc_nN_mi_y];
    fcs_nN_ave_mi_corrected{i} = fc_nN_mi_corrected;
end

figure(), hold on
for i = 1:length(fcs_nN_ave_mi_corrected)
    
    fc_x_mi = fcs_nN_ave_mi_corrected{i}(:,1);
    fc_y_mi = fcs_nN_ave_mi_corrected{i}(:,2);
    
    cp_idx = cp_idx_ave_mi(i);
    
    plot(fc_x_mi, fc_y_mi, 'k')
    scatter(fc_x_mi(cp_idx), fc_y_mi(cp_idx), 'g', 'x')
end
plot(fcs_nN_ave_mi_corrected{1}(:,1), fcs_nN_ave_mi_corrected{1}(:,2), 'r', 'LineWidth', 2)
scatter(fcs_nN_ave_mi_corrected{1}(cp_idx_ave_mi(1), 1), fcs_nN_ave_mi_corrected{1}(cp_idx_ave_mi(1), 2), 'r', 'x', 'LineWidth', 2)
xlabel('Tip-sample separation (nm)', 'FontSize', 13)
ylabel('Force (nN)', 'FontSize', 13)
title('Force curves (MI) - aligned on contact points', 'FontSize', 14)
pbaspect([2 1 1])
hold off

%% remove same area of indentation as for cp aligned fcs

fcs_nN_ave_mi_final = cell(size(fcs_nN_ave_mi_corrected));
fcs_pN_ave_mi_final = cell(size(fcs_nN_ave_mi_corrected));

cp_idx_ave_mi_final = zeros(size(cp_idx_ave_mi));

for i = 1:length(fcs_nN_ave_mi_corrected);

    fc = fcs_nN_ave_mi_corrected{i};
    
    cutoff_iv_nm = (-1) * cutoff_iv_nm_rb(i);
    fc_averaged_x    = fc(:,1);
    fc_averaged_nN_y = fc(:,2);
    fc_averaged_pN_y = fc_averaged_nN_y * 1000;
    
    idx = find(fc_averaged_x == cutoff_iv_nm);
    idx_baseline = find(fc_averaged_x == (fc_averaged_x(end)-20));
    
    if isempty(idx) == 1
        idx = 1;
    end
    fc_averaged_x_cropped = fc_averaged_x(idx:idx_baseline); 
    
    cpidx_final = find(fc_averaged_x_cropped == 0);
    cp_idx_ave_mi_final(i) = cpidx_final;
    
    fc_averaged_nN_y_cropped = fc_averaged_nN_y(idx:idx_baseline); 
    fc_averaged_pN_y_cropped = fc_averaged_pN_y(idx:idx_baseline); 
    
    
    fcs_nN_ave_mi_final{i} = [fc_averaged_x_cropped, fc_averaged_nN_y_cropped];
    fcs_pN_ave_mi_final{i} = [fc_averaged_x_cropped, fc_averaged_pN_y_cropped];
    
end

figure(), hold on
for i = 1:length(fcs_nN_ave_mi_final)
    
    fc_x_mi = fcs_nN_ave_mi_final{i}(:,1);
    fc_y_mi = fcs_nN_ave_mi_final{i}(:,2);
    
    cp_idx = cp_idx_ave_mi_final(i);
    
    plot(fc_x_mi, fc_y_mi, 'k')
    scatter(fc_x_mi(cp_idx), fc_y_mi(cp_idx), 'g', 'x')
end
plot(fcs_nN_ave_mi_final{1}(:,1), fcs_nN_ave_mi_final{1}(:,2), 'r', 'LineWidth', 2)
scatter(fcs_nN_ave_mi_final{1}(cp_idx_ave_mi_final(1), 1), fcs_nN_ave_mi_final{1}(cp_idx_ave_mi_final(1), 2), 'r', 'x', 'LineWidth', 2)
xlabel('Tip-sample separation (nm)', 'FontSize', 13)
ylabel('Force (nN)', 'FontSize', 13)
title('Force curves (MI) - removed indentation past threshold', 'FontSize', 14)
pbaspect([2 1 1])
hold off
    
%% Apply Hertz model to each force curve from the averaged contact point for fcs averaged by max indentation  

ym_mi_MPa = zeros(size(fcs_nN_ave_mi_final));

for i = 1:length(fcs_nN_ave_mi_final)
    
    fc_nN_mi = fcs_nN_ave_mi_final{i};
    cpidx    = cp_idx_ave_mi_final(i);  
    
    ym_MPa       = HertzModel(fc_nN_mi, cpidx, HertzIndentation_nm); 
    ym_mi_MPa(i) = ym_MPa;
    
end

%% True height and E_eff calculated from max indentation force curves

% make arrays for plotting

th_mirror    = fliplr(th_rb_ave(2:end));
ym_mi_mirror = fliplr(ym_mi_MPa(2:end));

th_plot    = [th_mirror, th_rb_ave];
ym_mi_plot = [ym_mi_mirror, ym_mi_MPa];

% plot rotationally averaged E_eff as calculated from: (1) appling Hertz 
% model to each force curve individually and then rotationally averaging 
% the results; and (2) aligning the force curves by maximum indentation
% (rather than contact point), averaging them, finding their contact point,
% and then applying the Hertz model. If the shape of the rotationally
% averaged plots are the same, we know that the contact point determination
% is not unduly affecting the results.
figure(), plot(rv_plot, ym_plot, 'Color', c1, 'LineWidth', 2)
hold on
plot(rv_plot, ym_mi_plot, 'k--', 'LineWidth', 2)
axis([-100 100 0 6])
xlabel('Radial distance (nm)', 'FontSize', 13)
ylabel('E_{eff} (MPa)', 'FontSize', 13);
legend('Hertz model applied to fcs individually', 'Hertz model applied to averaged fcs (MI)', 'Location', 'SouthOutside')
title(strcat('Rotationally averaged E_{eff} (n=', num2str(numbpores), ')'));
pbaspect([2 1 1])

figure(), plot(rv_plot, hd_plot, 'Color', c2, 'LineWidth', 2)
hold on
plot(rv_plot, th_plot, 'Color', c3, 'LineStyle', '--', 'LineWidth', 2)
axis([-100 100 0 60])
xlabel('Radial distance (nm)', 'FontSize', 13)
ylabel('Height (nm)', 'FontSize', 13);
set(gca, 'FontSize', 13)
title(strcat('Rotationally averaged height (n=', num2str(numbpores), ')'));
legend('Height at max indentation', 'True height', 'Location', 'SouthOutside')
pbaspect([2 1 1])

%% Now wish to repeat but for pores individually
% 1. Average the height data by pore

hd_pp_ave = cell(1, numbpores);
for i = 1:length(hd)
    for j = 1:length(hd{i})
        
        hd_pp_bin  = hd{i}{j};
        hd_pp_mean = mean(hd_pp_bin(:));
        
        hd_pp_ave{i}(j) = hd_pp_mean;
        
    end        
end


%% 2. Average the YM values by pore

ym_pp_ave         = cell(1, numbpores);
cc_pp_sum         = cell(1, numbpores);

for i = 1:length(ym_iv_ymidx)
    for j = 1:length(ym_iv_ymidx{i})
    
    cc_pp_bin_sum = sum(cc_iv_ymidx{i}{j}(:));
    cc_pp_sum{i}(j) = cc_pp_bin_sum;
        
    ym_pp_bin = ym_iv_ymidx{i}{j};
    ym_pp_bin_sum = sum(ym_pp_bin(:));
    ym_pp_ave{i}(j)  = ym_pp_bin_sum/cc_pp_bin_sum;
    
    end
    
end

%% 3. Average the force curves by pore - re-organise force curves by radial bin
% first, need to do as above, i.e., go through and remove nans. But this
% time, without re-ordering

fc_pp = cell(size(fc_iv_ymidx));

length_memory_fc = zeros(length(rv), 1);

for i = 1:numbpores % number of pores
    for j = 1:length(fc_iv_ymidx{i}) % number of bins       
                
        count = cc_iv_ymidx{i}{j};
        numb_fcs = sum(count(:)); % this is number of force curves in this bin of this pore
                
        % as cropped matrices can vary in size, so can the length of the
        % same radial bins between pores. Therefore must keep track of how
        % many pixels are in each radial bin to ensure we don't add in
        % extra zeros which will affect the averaging process.
        previous_length = length_memory_fc(j);
        length_memory_fc(j) = previous_length + numb_fcs;
        
        current_length = length_memory_fc(j);
        
        % as some force curves have been removed and replaced with an nan,
        % we need an fc_nan_counter to track them, and adjust the indexing
        % so that we only pull-out the real force curves
        fc_nan_counter = 0;
        
        if numb_fcs > 0 % if there are fcs in the bin
            
            for n = 1:length(count); % for length of pixels in bin (binned or not)
                
                if count(n) == 1 % and if the force curve exists           

                    fc_current = fc_iv_ymidx{i}{j}(n);
                    fc_pp{i}{j}(1, (n - fc_nan_counter)) = fc_current; % hd re-organised into radial bin positions
                
                else % if no force curve exists, the counter increases by 1, so that in the next iteration fc_current will be placed in idx n-1, thereby not skipping a position                 
                    fc_nan_counter = fc_nan_counter + 1;
                    continue
                end
                
            end
            
        end

    end
           
end

%% 3. (contd.) now average force curves

display('Averaging force curves by pore...')

% pre-allocate
fcs_pp_nN_ave = cell(1, numbpores);
fcs_pp_pN_ave = cell(1, numbpores);
ind_pp_rb_ave = cell(1, numbpores);
th_pp_ave     = cell(1, numbpores);

for i = 1:numbpores
    for j = 1:length(fc_pp{i})
        
        fcs_pp_to_ave = fc_pp{i}{j};
        [~, fc_pp_averaged, ~, ~] = NanoMechFuncs.Average_forcecurves(fcs_pp_to_ave, binwidth_nm, threshold);
        
        if isempty(fc_pp_averaged) == 0 % if not all force curves in radial bin removed (highly unlikely all will be removed)
        
            cutoff_iv_nm = (-1) * cutoff_iv_nm_rb(j);

            xdata    = fc_pp_averaged(:,1);
            ydata_nN = fc_pp_averaged(:,2);
            ydata_pN = ydata_nN .* 1000;

            idx          = find(xdata == cutoff_iv_nm);
            idx_baseline = find(fc_averaged_x == (fc_averaged_x(end)-20));
            if isempty(idx) == 1
                idx = 1;
            end


            xdata_cropped    = xdata(idx:idx_baseline);
            ydata_nN_cropped = ydata_nN(idx:idx_baseline);
            ydata_pN_cropped = ydata_pN(idx:idx_baseline);

            fcs_pp_nN_ave{i}{j} = [xdata_cropped, ydata_nN_cropped];
            fcs_pp_pN_ave{i}{j} = [xdata_cropped, ydata_pN_cropped];

            indentation_pp   = abs(min(xdata_cropped));
            ind_pp_rb_ave{i}(j) = indentation_pp;

            trueheight_pp   = hd_pp_ave{i}(j) + indentation_pp;
            th_pp_ave{i}(j) = trueheight_pp;
            
        end
        
    end
    
end

%% 4. Mirror the hd, ym, and th data for plotting

hd_pp_plot = cell(1, numbpores);
ym_pp_plot = cell(1, numbpores);
th_pp_plot = cell(1, numbpores);

for i = 1:numbpores
    
    hd_pp_ave_i = hd_pp_ave{i};
    ym_pp_ave_i = ym_pp_ave{i};
    th_pp_ave_i = th_pp_ave{i};
    
    hd_pp_mirror = fliplr(hd_pp_ave_i(2:end));
    ym_pp_mirror = fliplr(ym_pp_ave_i(2:end));
    th_pp_mirror = fliplr(th_pp_ave_i(2:end));

    
    hd_pp_plot{i} = [hd_pp_mirror, hd_pp_ave_i];
    ym_pp_plot{i} = [ym_pp_mirror, ym_pp_ave_i];
    th_pp_plot{i} = [th_pp_mirror, th_pp_ave_i];

end

%% Plot E_eff with scatter from all pores

grey = 0.6;

figure(),
hold on
for i=1:length(ym_pp_plot)

        ym_data = ym_pp_plot{i}(:);
        scatter(rv_plot, ym_data, 'MarkerEdgeColor', grey*[1 1 1], 'LineWidth', 0.5);  
    
end
hold on
errorbar(rv_plot, ym_plot, ym_rb_halfstd_plot, 'Color', c1, 'LineWidth', 2)
xlabel('Radial position (nm)', 'FontSize', 13)
ylabel('E_{eff} (MPa)', 'FontSize', 13)
title(['E_{eff} (MPa)', ' (n=', num2str(numbpores), ')'], 'FontSize', 14)
xlim([-100 100])
ylim([0 6])
pbaspect([2 1 1])

%% Plot all heights (black) with true heights (blue) and averages (red and green respectively)

figure(),
hold on
for i = 1:length(hd_pp_plot)
    plot(rv_plot, hd_pp_plot{i}, 'k', 'LineWidth', 1)
    hold on
    plot(rv_plot, th_pp_plot{i}, 'b', 'LineWidth', 1)
end
hold on
plot(rv_plot, hd_plot, 'r', 'LineWidth', 2)
hold on
plot(rv_plot, th_plot, 'g', 'LineWidth', 2)
xlabel('Radial position (nm)', 'FontSize', 13)
ylabel('Height (nm)', 'FontSize', 13)
title('Height at MI (black; red=ave); true height (blue; green=ave)')
xlim([-100 100])
ylim([0 60])
pbaspect([2 1 1])


%% calculate the percentage of force curves removed by bin and in total and plot height and ym to check

cc_lengths_all  = sum(cc_rb_length(:));
cc_sum_all      = sum(cc_rb_sum(:));
percent_all_bin = 100 - ((cc_sum_all/cc_lengths_all)*100);

display(percent_rb_bin)
display(percent_all_bin)

%% Stiffness Curves - overall average

% This is for finding the stiffness curves, and should always be set to
% 'yes'.
Invert_Matlab_Derivative = 'Yes';
% with force curves aligned by contact point
[scs_pN_ave, ~] = NanoMechFuncs.DerivativeForceCurve_2(fcs_pN_ave, Invert_Matlab_Derivative);
% with force curves aligned by maximum indentation
[scs_pN_ave_mi, ~] = NanoMechFuncs.DerivativeForceCurve_2(fcs_pN_ave_mi_final, Invert_Matlab_Derivative);

%% Plot force curves aligned by contact point or max indentation

figure(),
hold on
for i = 1:length(fcs_pN_ave)
    
    fc_x = fcs_pN_ave{i}(:,1);
    fc_y = fcs_pN_ave{i}(:,2);
    
    plot(fc_x, fc_y, 'k', 'LineWidth', 1);
    
end
plot(fcs_pN_ave{1}(:,1), fcs_pN_ave{1}(:,2), 'r', 'LineWidth', 2)
title('Force curves (aligned by CP)', 'FontSize', 15)
xlabel('Tip-sample separation (nm)', 'FontSize', 13)
ylabel('Force (pN)', 'FontSize', 13)
pbaspect([2 1 1])
hold off

figure(), hold on
for i = 1:length(fcs_pN_ave_mi_final)
    
    plot(fcs_pN_ave_mi_final{i}(:,1), fcs_pN_ave_mi_final{i}(:,2), 'k', 'LineWidth', 1)
    
end
plot(fcs_pN_ave_mi_final{1}(:,1), fcs_pN_ave_mi_final{1}(:,2), 'r', 'LineWidth', 2)
title('Force curves (aligned by MI)', 'FontSize', 15)
xlabel('Tip-sample separation (nm)', 'FontSize', 13)
ylabel('Force (pN)', 'FontSize', 13)
pbaspect([2 1 1])
hold off

%% Plot stiffness curves aligned by contact point or max indentation

figure(),
hold on
for i = 1:length(scs_pN_ave)
    
    sc_x = scs_pN_ave{i}(:,1);
    sc_y = scs_pN_ave{i}(:,2);
    
    plot(sc_x, sc_y, 'k', 'LineWidth', 1);
    
end
plot(scs_pN_ave{1}(:,1), scs_pN_ave{1}(:,2), 'r', 'LineWidth', 2)
title('Stiffness (aligned by CP)', 'FontSize', 15)
xlabel('Tip-sample separation (nm)', 'FontSize', 13)
ylabel('Stiffness (pN/nm)', 'FontSize', 13)
pbaspect([2 1 1])
hold off

figure(),
hold on
for i = 1:length(scs_pN_ave_mi)
    
    sc_x_mi = scs_pN_ave_mi{i}(:,1);
    sc_y_mi = scs_pN_ave_mi{i}(:,2);
    
    plot(sc_x_mi, sc_y_mi, 'k', 'LineWidth', 1);
    
end
plot(scs_pN_ave_mi{1}(:,1), scs_pN_ave_mi{1}(:,2), 'r', 'LineWidth', 2)
title('Stiffness (aligned by MI)', 'FontSize', 15)
xlabel('Tip-sample separation (nm)', 'FontSize', 13)
ylabel('Stiffness (pN/nm)', 'FontSize', 13)
pbaspect([2 1 1])
hold off

%% Heat map for fcs aligned and averaged by maximum indentation

% for plotting, set the edge of the true height to 0nm
th_value    = th_plot(end);
th_plot_0nm = th_plot - th_value;
hd_plot_0nm = hd_plot - th_value;

hd_rb_ave_adj = hd_rb_ave - th_value;

Height_Axis_Max = round(max(rv_plot));
Height_Axis_Min = round(min(rv_plot));

Height_Profile = sort([Height_Axis_Min: binwidth_nm :Height_Axis_Max]', 'ascend');

% create matrix with pN/nm values aligned radially and on averaged height
% at maximum indentation
scs_pN_ave_mi_adj = cell(1,length(scs_pN_ave_mi));
stiff_z_min = zeros(1,length(scs_pN_ave_mi));
    for i=1:length(scs_pN_ave_mi);
        
        stiff_z_min(i) = min(scs_pN_ave_mi{i}(:,1));
       
        scs_pN_ave_mi_adj{i}(:,1) = scs_pN_ave_mi{i}(:,1) + (hd_rb_ave_adj(i) - stiff_z_min(i));
        scs_pN_ave_mi_adj{i}(:,2) = scs_pN_ave_mi{i}(:,2);
        
    end
    
% Mirror the force curves for the plotting
scs_pN_ave_mi_adj_mirrored = fliplr(scs_pN_ave_mi_adj(2:end));
scs_pN_ave_mi_adj_hmap = horzcat(scs_pN_ave_mi_adj_mirrored, scs_pN_ave_mi_adj);

% create the stiffness heat map by finding the correct bin for each value
% of the stiffness curve.
stiff_pN_hmap_mi = zeros(length(Height_Profile),length(rv_plot));
for i=1:length(scs_pN_ave_mi_adj_hmap)
    for j=1:length(scs_pN_ave_mi_adj_hmap{i})
        Tip_Sample_Sep_Value        = scs_pN_ave_mi_adj_hmap{i}(j,1);
        Closest_Height              = abs(Height_Profile - Tip_Sample_Sep_Value);
        [rr_hm,~]                   = find(Closest_Height == min(Closest_Height), 1);
        stiff_pN_hmap_mi(rr_hm,i)   = scs_pN_ave_mi_adj_hmap{i}(j,2);
    end
end


% plot stiffness heat map in pN/nm, with a small interpolation to aid the
% plotting.
stiffhmap_pN_interp_mi = interp2((stiff_pN_hmap_mi), 5); % change from nN to pN

% My custom colour map 'Blue_Yellow_Red' is saved to the Matlab structure
% 'MyColourmaps'. First must load the structure, and then assign a variable
% to the desired colour map saved in the structure.
Structure     = load('MyColourMaps');
BYR_Colourmap = Structure.MyColourMaps.BYR;

% plot stiffness heatmap
max_intensity = max(stiffhmap_pN_interp_mi(:));
cHigh = max_intensity * 0.95;
cLow = 0;
figure(), imagesc(stiffhmap_pN_interp_mi, 'XData', [min(rv_plot) max(rv_plot)], 'YData', [min(Height_Profile) max(Height_Profile)])
set(gca, 'FontSize', 13)
set(gca, 'CLim', [cLow, cHigh]);
colormap(BYR_Colourmap); 
g = colorbar;
ylabel(g,'Stiffness (pN/nm)', 'FontSize', 13)
set(gca,'ydir','normal');
hold on
xlabel('Radial distance (nm)', 'FontSize', 13)
ylabel('Tip-Sample Separation (nm)', 'FontSize', 13)
title('Heat Map (aligned by MI) (pN/nm)')
hold on
plot(rv_plot, th_plot_0nm, '--','Color', [0 0 0], 'LineWidth', 3)
hold on
h = area(rv_plot, hd_plot_0nm, -80); % the -80 specifies the baseline value
h.FaceColor = [0.8 0.8 0.8];
h.EdgeColor = [0.8 0.8 0.8];
axis([-100 100 -50 50])
pbaspect([2 1 1])

%% Heat maps - aligned by contact point

% for plotting, set the edge of the true height to 0nm
th_value    = th_plot(end);
th_plot_0nm = th_plot - th_value;
hd_plot_0nm = hd_plot - th_value;

hd_rb_ave_adj = hd_rb_ave - th_value;

Height_Axis_Max = round(max(rv_plot));
Height_Axis_Min = round(min(rv_plot));
Height_Profile  = sort([Height_Axis_Min: binwidth_nm :Height_Axis_Max]', 'ascend');

% create matrix with pN/nm values aligned radially and on averaged height
% at maximum indentation
scs_pN_ave_adj = cell(1,length(scs_pN_ave));
stiff_z_min = zeros(1,length(scs_pN_ave));
    for i=1:length(scs_pN_ave);
        
        stiff_z_min(i) = min(scs_pN_ave{i}(:,1));
       
        scs_pN_ave_adj{i}(:,1) = scs_pN_ave{i}(:,1) + (hd_rb_ave_adj(i) - stiff_z_min(i));
        scs_pN_ave_adj{i}(:,2) = scs_pN_ave{i}(:,2);
    end
    
% Mirror the force curves for the plotting
scs_pN_ave_adj_mirrored = fliplr(scs_pN_ave_adj(2:end));
scs_pN_ave_adj_hmap = horzcat(scs_pN_ave_adj_mirrored, scs_pN_ave_adj);

% create the stiffness heat map by finding the correct bin for each value
% of the stiffness curve.
stiff_pN_hmap = zeros(length(Height_Profile),length(rv_plot));
for i=1:length(scs_pN_ave_adj_hmap)
    for j=1:length(scs_pN_ave_adj_hmap{i})
        Tip_Sample_Sep_Value        = scs_pN_ave_adj_hmap{i}(j,1);
        Closest_Height              = abs(Height_Profile - Tip_Sample_Sep_Value);
        [rr_hm,~]                   = find(Closest_Height == min(Closest_Height), 1);
        stiff_pN_hmap(rr_hm,i) = scs_pN_ave_adj_hmap{i}(j,2);
    end
end

% plot stiffness heat map in pN/nm, with a small interpolation to aid the
% plotting.
stiffhmap_pN_interp = interp2((stiff_pN_hmap), 5); % change from nN to pN

% My custom colour map 'Blue_Yellow_Red' is saved to the Matlab structure
% 'MyColourmaps'. First must load the structure, and then assign a variable
% to the desired colour map saved in the structure.
Structure = load('MyColourMaps');
BYR_Colourmap = Structure.MyColourMaps.BYR;

% plot stiffness heatmap
max_intensity = max(stiffhmap_pN_interp(:));
cHigh = max_intensity * 0.95;
cLow = 0;
figure(), imagesc(stiffhmap_pN_interp, 'XData', [min(rv_plot) max(rv_plot)], 'YData', [min(Height_Profile) max(Height_Profile)])
set(gca, 'FontSize', 13)
set(gca, 'CLim', [cLow, cHigh]);
colormap(BYR_Colourmap); 
g = colorbar;
ylabel(g,'Stiffness (pN/nm)', 'FontSize', 13)
set(gca,'ydir','normal');
hold on
xlabel('Radial distance (nm)', 'FontSize', 13)
ylabel('Tip-Sample Separation (nm)', 'FontSize', 13)
title('Heat Map (aligned by CP) (pN/nm)')
hold on
plot(rv_plot, th_plot_0nm, '--','Color', [0 0 0], 'LineWidth', 3)
hold on
h = area(rv_plot, hd_plot_0nm, -80); % the -80 specifies the baseline value
h.FaceColor = [0.8 0.8 0.8];
h.EdgeColor = [0.8 0.8 0.8];
axis([-100 100 -50 50])
pbaspect([2 1 1])

%% Now, repeat for pores individually - Stiffness Curves: per pore average

scs_pp_pN_ave               = cell(size(fcs_pp_pN_ave));
stiff_pp_pN_hmap_cell       = cell(1, numbpores);
stiffhmap_pp_pN_interp_cell = cell(1, numbpores);
th_pp_plot_0nm_cell         = cell(1, numbpores);
hd_pp_plot_0nm_cell         = cell(1, numbpores);

Invert_Matlab_Derivative = 'Yes';

for i = 1:length(fcs_pp_pN_ave)
    
    fcs_pp_pN_ave_cell = fcs_pp_pN_ave{i};
        
    [scs_pp_pN_ave_cell, ~] = DerivativeForceCurve_2(fcs_pp_pN_ave_cell, Invert_Matlab_Derivative);
    scs_pp_pN_ave{i} = scs_pp_pN_ave_cell;
    
end

%% Make heatmaps (CP) for each pore individually and save

for n = 1:numbpores

    th_pp_plot_array = th_pp_plot{n};
    hd_pp_plot_array = hd_pp_plot{n};
    
    hd_pp_rb_ave = hd_pp_ave{n};
    
    % for plotting, set the edge of the true height to 0nm
    th_pp_value    = th_pp_plot_array(end);
    th_pp_plot_0nm = th_pp_plot_array - th_pp_value;
    hd_pp_plot_0nm = hd_pp_plot_array - th_pp_value;
    
    % save these out for plotting in final script
    th_pp_plot_0nm_cell{n} = th_pp_plot_0nm;
    hd_pp_plot_0nm_cell{n} = hd_pp_plot_0nm;

    hd_pp_rb_ave_adj = hd_pp_rb_ave - th_pp_value;

    Height_Axis_Max = round(max(rv_plot));
    Height_Axis_Min = round(min(rv_plot));

    Height_Profile = sort([Height_Axis_Min: binwidth_nm :Height_Axis_Max]', 'ascend');

    scs_pp_pN_ave_cell_array = scs_pp_pN_ave{n};
    scs_pp_pN_ave_adj = cell(1,length(scs_pp_pN_ave_cell_array));

    stiff_z_pp_min = zeros(1,length(scs_pp_pN_ave_cell_array));
        for i=1:length(scs_pp_pN_ave_cell_array);
            
            if isempty(scs_pp_pN_ave_cell_array{i}) == 0

                stiff_z_pp_min(i) = min(scs_pp_pN_ave_cell_array{i}(:,1));

                scs_pp_pN_ave_adj{i}(:,1) = scs_pp_pN_ave_cell_array{i}(:,1) + (hd_pp_rb_ave_adj(i) - stiff_z_pp_min(i));
                scs_pp_pN_ave_adj{i}(:,2) = scs_pp_pN_ave_cell_array{i}(:,2);
                
            end
        end

    % Mirror the force curves for the plotting
    scs_pp_pN_ave_adj_mirrored = fliplr(scs_pp_pN_ave_adj(2:end));
    scs_pp_pN_ave_adj_hmap = horzcat(scs_pp_pN_ave_adj_mirrored, scs_pp_pN_ave_adj);

    % create the stiffness heat map by finding the correct bin for each value
    % of the stiffness curve.
    stiff_pp_pN_hmap = zeros(length(Height_Profile),length(rv_plot));
    for i=1:length(scs_pp_pN_ave_adj_hmap)
        for j=1:length(scs_pp_pN_ave_adj_hmap{i})
            Tip_Sample_Sep_Value        = scs_pp_pN_ave_adj_hmap{i}(j,1);
            Closest_Height              = abs(Height_Profile - Tip_Sample_Sep_Value);
            [rr_hm,~]                   = find(Closest_Height == min(Closest_Height), 1);
            stiff_pp_pN_hmap(rr_hm,i) = scs_pp_pN_ave_adj_hmap{i}(j,2);
        end
    end
    
    % save out for plotting in final script
    stiff_pp_pN_hmap_cell{n} = stiff_pp_pN_hmap;

    % plot stiffness heat map in pN/nm, with a small interpolation to aid the
    % plotting.
    stiffhmap_pp_pN_interp         = interp2((stiff_pp_pN_hmap), 5); 
    % save these out for plotting in final script
    stiffhmap_pp_pN_interp_cell{n} = stiffhmap_pp_pN_interp;

end
    
%% Store data in a structure and then save for plotting in the final script

if save_out_data_structure == 1

    display('Saving data...')

    % create a writable mat file
    FullFileOutput = fullfile(OutputFolder, strcat(GenericFileName, ' - NanomechanicalProcessedData'));
    NanomechanicalProcessedData = matfile(strcat(FullFileOutput, '.mat'), 'Writable', true);

    % save data out for the pores individually
    PP_Processed.fcs_pp_nN_ave                        = fcs_pp_nN_ave;
    PP_Processed.fcs_pp_pN_ave                        = fcs_pp_pN_ave;
    PP_Processed.scs_pp_pN_ave                        = scs_pp_pN_ave;
    PP_Processed.hd_pp_ave                            = hd_pp_ave;
    PP_Processed.ym_pp_ave                            = ym_pp_ave;
    PP_Processed.th_pp_ave                            = th_pp_ave;
    PP_Processed.hd_pp_plot                           = hd_pp_plot;
    PP_Processed.ym_pp_plot                           = ym_pp_plot;
    PP_Processed.th_pp_plot                           = th_pp_plot;
    PP_Processed.rv_plot                              = rv_plot;
    % put in separate field the variables required for plotting the heat maps
    PP_Processed.hmapplot.stiff_pp_pN_hmap_cell       = stiff_pp_pN_hmap_cell;
    PP_Processed.hmapplot.stiffhmap_pp_pN_interp_cell = stiffhmap_pp_pN_interp_cell;
    PP_Processed.hmapplot.th_pp_plot_0nm_cell         = th_pp_plot_0nm_cell;
    PP_Processed.hmapplot.hd_pp_plot_0nm_cell         = hd_pp_plot_0nm_cell;

    % write to the matfile
    NanomechanicalProcessedData.PP_Processed = PP_Processed;
    
    % save data for all averaged pores
    RB_Processed.fcs_nN_ave                      = fcs_nN_ave;
    RB_Processed.fcs_pN_ave                      = fcs_pN_ave;
    RB_Processed.scs_pN_ave                      = scs_pN_ave;
    RB_Processed.fcs_pN_ave_mi                   = fcs_pN_ave_mi_final ;
    RB_Processed.scs_pN_ave_mi                   = scs_pN_ave_mi;
    RB_Processed.hd_rb_ave                       = hd_rb_ave;
    RB_Processed.ym_rb_ave                       = ym_rb_ave;
    RB_Processed.th_rb_ave                       = th_rb_ave;
    RB_Processed.hd_plot                         = hd_plot;
    RB_Processed.ym_plot                         = ym_plot;
    RB_Processed.ym_mi_plot                      = ym_mi_plot;
    RB_Processed.th_plot                         = th_plot;
    RB_Processed.rv_plot                         = rv_plot;
    % save out standard deviations
    RB_Processed.ym_rb_std_plot                  = ym_rb_std_plot;
    RB_Processed.ym_rb_halfstd_plot              = ym_rb_halfstd_plot;
    RB_Processed.hd_rb_std_plot                  = hd_rb_std_plot;
    RB_Processed.hd_rb_halfstd_plot              = hd_rb_halfstd_plot;
    % and percentage of binning
    RB_Processed.percent_rb_bin                  = percent_rb_bin;
    RB_Processed.percent_all_bin                 = percent_all_bin;
    % put in separate field the variables required for plotting the heat maps
    RB_Processed.hmapplot.stiff_pN_hmap          = stiff_pN_hmap;
    RB_Processed.hmapplot.stiffhmap_pN_interp    = stiffhmap_pN_interp;
    RB_Processed.hmapplot.stiff_pN_hmap_mi       = stiff_pN_hmap_mi;
    RB_Processed.hmapplot.stiffhmap_pN_interp_mi = stiffhmap_pN_interp_mi;
    RB_Processed.hmapplot.th_plot_0nm            = th_plot_0nm;
    RB_Processed.hmapplot.hd_plot_0nm            = hd_plot_0nm;
    
    
    % write to the matfile
    NanomechanicalProcessedData.RB_Processed = RB_Processed;
    % carry through the cropped image and YM data
    % write to the matfile
    NanomechanicalProcessedData.radiallybinned = radiallybinned;
    NanomechanicalProcessedData.entirematrices = entirematrices;
    NanomechanicalProcessedData.matrices       = matrices;

end


