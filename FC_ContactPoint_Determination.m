% ContactPoint_12 simply crops the force curve from 1 interval beyond the
% unique transition point, to the end. If this means the entire force curve
% is taken for the piecewise function, then this is allowed.

% However, rather than the reference being Factor*min(StDv), we now take
% the median of the StDv, rather than the min.

function [Force_Curve, Force_Curve_Corrected, Contact_Index, Indentation] = FC_ContactPoint_Determination(Force_Curve, Factor, Smoothe, Number_moving_ave_filter)

% Check the gradient from the first point to the last. If the force curve
% has been input in the wrong format, it is switched around.

x1 = Force_Curve(1,1);
x2 = Force_Curve(end,1);

% if the force curve in the arrays is given in the opposite direction to
% those from QNM files, it will be upended.

if x1 > x2
    Force_Curve = flipud(Force_Curve);
end

% Rename xData and yData, as will know Force Curve is in correct direction.
xTrace_Unfiltered = Force_Curve(:,1);
yTrace_Unfiltered = Force_Curve(:,2);

%% Smoothe data

% if decide to use a moving average filter
if length(Smoothe) == 3
    xTrace = Force_Curve(:,1);
    yTrace = smooth(Force_Curve(:,2),Number_moving_ave_filter);
else
    xTrace = Force_Curve(:,1);
    yTrace = Force_Curve(:,2);
end



Force_Curve = [xTrace,yTrace];
Data_Length = length(Force_Curve);


%% Pre-allocate cell arrays
Count                 = zeros(1,18);
Interval_Find         = cell(1,18);
Interval              = cell(1,18);
Mean                  = cell(1,length(Interval));
StDv                  = cell(1,length(Interval));
Mean_Diff             = cell(1,length(Interval));
True_False_Count      = cell(1,length(Interval));
True_False_Diff       = cell(1,length(Interval));
Force_Curve_Corrected = zeros(Data_Length,2);

%%%%%% Mean and StDv of intervals
Interval_Length = length(Interval); % 18

xTrace_min = min(xTrace);
xTrace_max = max(xTrace);
xTrace_d   = abs(xTrace_max - xTrace_min);

Interval_Find{1} = zeros(1,21);
Interval{1}      = zeros(1,21);
for n=1:length(Interval_Find);
    Interval_Find{n}(1,1) = xTrace_min; % ensure first element is always xTrace_min
    Interval{n}(1,1)      = 1; % ensure first element is always 1
    Int_d                 = (xTrace_d/(length(Interval_Find{n})-1));
    for j=1:length(Interval_Find{n})-1;
        Interval_Find{n}(1,j+1) = (Int_d * j) + xTrace_min;
        xTrace_Closest          = abs(xTrace - Interval_Find{n}(j+1));
        Interval{n}(j+1)        = find(xTrace_Closest == min(xTrace_Closest), 1 );
    end
    if n == Interval_Length;
        break
    end
    Interval_Find{n+1} = zeros(1,length(Interval_Find{n})-1);
    Interval{n+1}      = zeros(1,length(Interval{n})-1);
end

%% Calculate the mean and standard deviation for each interval and store
% these values in their corresponding cell arrays.
for n=1:length(Interval)
    Interval_Data = cell(1,length(Interval{n})-1);
    Mean{n} = zeros(1,length(Interval{n})-1);
    StDv{n} = zeros(1,length(Interval{n})-1);
    for j=1:length(Interval{n})-1
        Interval_Data{j} = Force_Curve((Interval{n}(j):(Interval{n}(j+1))), 2);
        Mean{n}(j) = mean(Interval_Data{j});
        StDv{n}(j) =  std(Interval_Data{j});
    end
end

%% Absolute Difference Between the Mean Values
for n=1:length(Interval)
    Mean_Diff{n} = zeros(1,length(Mean{n})-1);
    for j=1:length(Mean{n})-1
        Mean_Diff{n}(j) = abs(Mean{n}(j)-Mean{n}(j+1));
    end
end

% Find intervals with a unique transition point
Reference = zeros(1,length(Interval));
for n=1:length(Interval)
    Reference(n) = median(StDv{n})*Factor;
    True_False_Count{n} = zeros(1,length(Mean_Diff{n}));
    True_False_Diff{n}  = zeros(1,length(True_False_Count{n})-1);
    for j=1:length(Mean_Diff{n})
        if Mean_Diff{n}(j) <= Reference(n)
            True_False_Count{n}(j) = 0;
        else
            True_False_Count{n}(j) = 1;
        end
    end
        for m=1:length(True_False_Count{n})-1
            True_False_Diff{n}(m) = abs(True_False_Count{n}(m) - True_False_Count{n}(m+1));
        end
        Count(n) = sum(True_False_Diff{n});
end



%%

% The variable 'Iteration' is the number of iterations required to get 1
% unique transition point. This is stored and later used to crop the force
% curve for the fitting of the Piecewise function.

% It is assumed that the contact point will be found in the contact-region
% half of the force curve (i.e. from 0 to 0.5), and if it is not, the while loop continues, but
% with each iteration it adds more of the contact region for the fitting of the Piecewise
% function.

% for position as a fraction of cropped force curve length (cropped for
% fitting of piecewise function), use:
% Relative_Position_Crop   = 1;
% for position as a fraction of entire FC length, use:
Relative_Position = 1;
Attempt = 1;
while Relative_Position >= 0.5

    % If no unique transition point was found, the contact point is given as
    % NaN and the while loop breaks.
        if length(find(Count~=1)) == length(Count)
        Contact_Index                 = NaN;
        Contact_Points(1,1)           = NaN;
        Contact_Points(1,2)           = NaN;
        Contact_Points_Corrected(1,1) = NaN;
        Contact_Points_Corrected(1,2) = NaN;
        Force_Curve                   = NaN(Data_Length,2);
        Force_Curve_Corrected         = NaN(Data_Length,2);
        Baseline                      = NaN; 
        Indentation                   = NaN;
            break
        end

    [~,Index_Count] = find(Count == 1);

    % If while loop attempts to find the contact point again, but there are no
    % more possible iterations, the loop breaks.
        if Attempt > length(Index_Count)
        Contact_Index                 = NaN;
        Contact_Points(1,1)           = NaN;
        Contact_Points(1,2)           = NaN;
        Contact_Points_Corrected(1,1) = NaN;
        Contact_Points_Corrected(1,2) = NaN;
        Force_Curve                   = NaN(Data_Length,2);
        Force_Curve_Corrected         = NaN(Data_Length,2);
        Baseline                      = NaN; 
        Indentation                   = NaN;
            break
        end

    Iteration   = Index_Count(Attempt);
    Intervals   = Interval{Iteration};

    % If the iterations ran until the force curve is split into thirds, and
    % still there is not one unique transition point, the pixel is removed.
        if Iteration == length(Interval)
        Contact_Index                 = NaN;
        Contact_Points(1,1)           = NaN;
        Contact_Points(1,2)           = NaN;
        Contact_Points_Corrected(1,1) = NaN;
        Contact_Points_Corrected(1,2) = NaN;
        Force_Curve                   = NaN(Data_Length,2);
        Force_Curve_Corrected         = NaN(Data_Length,2);
        Baseline                      = NaN; 
        Indentation                   = NaN;
            break
        end

    % Find Transition Point. First, find the position (column) of all the
    % ones in the True_False_Count array. There will be more than one answer,
    % as this array contains number-of-intervals - 1 elements. Therefore, the second step, is
    % to take the maximum index value, which corresponds to the one closest to
    % the transition point.
    [~, Index_c]            = find(True_False_Count{Iteration} == max(True_False_Count{Iteration}));
    Transition_Point_Crop   = max(Index_c);
    Unique_Transition_Point = Transition_Point_Crop + 1;



    % Cropping_Point          = Transition_Point_Crop + 3;

    % If the Transition Point Crop is less than 1, the pixel is removed. This is 
    % because the next step is to crop the force curve from this point for the 
    % fitting of the Piecewise function. If the value is less than 1, you will
    % exceed the Matrix Dimensions. A Transition_Point_Crop of 1 is allowed,
    % although this means the fitting is done to the entire force curve.
    % Similarly, if the UTP is found to be a value greater than the
    % number of intervals, the pixel is removed.
        if Transition_Point_Crop < 1 || Unique_Transition_Point > length(Intervals)
        Contact_Index                 = NaN;
        Contact_Points(1,1)           = NaN;
        Contact_Points(1,2)           = NaN;
        Contact_Points_Corrected(1,1) = NaN;
        Contact_Points_Corrected(1,2) = NaN;
        Force_Curve                   = NaN(Data_Length,2);
        Force_Curve_Corrected         = NaN(Data_Length,2);
        Baseline                      = NaN;  
        Zmin                          = NaN;
            break
        end



    xdata     = Force_Curve(Intervals(Transition_Point_Crop):Intervals(end), 1);
    ydata     = Force_Curve(Intervals(Transition_Point_Crop):Intervals(end), 2);




    % Shift the cropped data so it begins from zero. This aids the fitting
    % function.
    xdata_one     = xdata(1);
    xdata_shifted = xdata - xdata_one;

    %%%==================== Boundary Conditions ============================%%%
    % Take the gradient of a line from the first point to the last point of the
    % original force curve. This is then used as a boundary condition for the
    % fitting of the Piecewise function.

    % This block of code assigns the indices for the cropping point and the
    % unique transition point in the cropped data. This means the
    % Cropped_TPC_Index should always be 1, and the Cropped_UTP_Index will
    % usually be one interval increased (in index values).
    %TPC_Index                       = Intervals(Transition_Point_Crop);

    TPC                       = Intervals(Transition_Point_Crop);
    Cropped_Interval_Indices  = Intervals(Transition_Point_Crop:end);
    UTP                       = Intervals(Unique_Transition_Point);
    TPCE                      = Intervals(end);

    min_Cropped_Interval_Indices    = min(Cropped_Interval_Indices) - 1;
    Cropped_UTP_Index               = UTP - min_Cropped_Interval_Indices;
    Cropped_TPC_Index               = TPC - min_Cropped_Interval_Indices;
    Cropped_TPCE_Index              = TPCE - min_Cropped_Interval_Indices;

    Max_contact_region       = max(ydata(Cropped_TPC_Index:Cropped_UTP_Index));
    Min_contact_region       = min(ydata(Cropped_TPC_Index:Cropped_UTP_Index));
    Mean_contact_region      = mean(ydata(Cropped_TPC_Index:Cropped_UTP_Index));
    Mean_non_contact_region  = mean(ydata(Cropped_UTP_Index:end));

        if Cropped_TPCE_Index > length(ydata)   
        Contact_Index                 = NaN;
        Contact_Points(1,1)           = NaN;
        Contact_Points(1,2)           = NaN;
        Contact_Points_Corrected(1,1) = NaN;
        Contact_Points_Corrected(1,2) = NaN;
        Force_Curve                   = NaN(Data_Length,2);
        Force_Curve_Corrected         = NaN(Data_Length,2);
        Baseline                      = NaN;  
        Indentation                   = NaN;
            break
        end

    Mean_non_contact_cropped = mean(ydata(Cropped_UTP_Index:Cropped_TPCE_Index));

    % If there is a large dip (similar to a jump to contact) just before the
    % contact point, it is possible that this will be chosen as the UTP. To
    % check for this, here, we calculate the trend in the gradients of the
    % contact, and non-conatact regions of the cropped force curve. contact_m
    % is the slope from the intersection with the UTP to the end of the contact
    % region (but rather than using this force value, the mean force value of
    % the interval is taken). noncontact_m is the same but of the non-contact
    % region. If contact_m > noncontact_m (less negative than), then the while
    % loop goes to the next iteration.
    contact_m    = ((Mean_contact_region - ydata(Cropped_UTP_Index))/(xdata(1) - xdata(Cropped_UTP_Index)));
    noncontact_m = ((ydata(Cropped_UTP_Index) - Mean_non_contact_cropped)/(xdata(Cropped_UTP_Index) - xdata(end)));

    if contact_m > noncontact_m || contact_m > 0
        Attempt = Attempt + 1;
        continue
    end

    ub_y1 = Mean_contact_region;
    ub_y2 = Mean_non_contact_region;
    ub_x1 = xdata_shifted(Cropped_UTP_Index);
    ub_x2 = xdata_shifted(end);

    % ub_m. Point 1 is where the mean contact region force value meets the unique
    % transition point x value. Point 2 is where the mean non-contact force
    % value meets the end of the non-contact region's x value. ub_m is the
    % gradient between Points 1 and 2, and is denoted the least negative value
    % (i.e. the mildest) for the slope. As the slope is always negative, it
    % represents the upper boundary for the least squares regerssion fitting.
    ub_m = (ub_y1-ub_y2)/(ub_x1-ub_x2);
    ub_slope     = ub_m;

    lb_y1 = Max_contact_region;
    lb_y2 = Mean_non_contact_region;
    lb_x1 = xdata_shifted(Cropped_TPC_Index);
    lb_x2 = xdata_shifted(Cropped_UTP_Index); % x at UTP

    % lb_m. Point 1 is where the max force from the contact region meets the
    % first x value (which is 0; this point will almost always be the fist
    % value in the array). Point 2 is the mean non-contact force meeting the
    % UTP. lb_m is this slope, and it is later multiplied by 3.
    lb_m = (lb_y1-lb_y2)/(lb_x1-lb_x2);
    lb_slope     = lb_m*3;
    lb_slope_10  = lb_m*10;

    % baseline must be between the minimum value and 3 fifths up the cropped
    % force curve.
    lb_baseline  = min(ydata);
    ub_base      = abs((max(ydata) - min(ydata))/5);
    ub_baseline  = min(ydata) + (3*ub_base); % upper boundary for baseline is 3 fifths up the cropped force curve

    % to obtain the boundaries for the intersects, two imaginary lines have
    % been drawn. The line giving the lower boundary is from where the minimum
    % y of the contact region hits the UTP, to x = 0 at a gradient of ub_m (the
    % mildest gradient calculated for the slope previously). For the upper
    % boundary of the intersect, a line from where the max y value meets the
    % UTP, to x = 0, at a gradient of 10 times the lb_m (the steepest slope
    % calculated previously).
    lb_intersect = Min_contact_region - (ub_m * lb_x2);
    ub_intersect = Min_contact_region - (lb_slope_10 * lb_x2);   

    if lb_baseline > ub_baseline || lb_intersect > ub_intersect|| lb_slope > ub_slope
        Attempt = Attempt + 1;
        continue
    end
    %%%=================== End Boundary Conditions =========================%%%

    %%%===================== Initial Conditions ============================%%%

    ic_m        = lb_m;
    ic_c        = lb_y1 - (lb_m*lb_x1);
    ic_baseline = Mean_non_contact_region;

    %%%==================== End Initial Conditions =========================%%%

    % suppress the 'Optimization complete' message
    opts = optimset('Display', 'off');
    
    % Fit the Piecewise function.
    F = @(B,xdata_shifted) max(B(1), B(2)+ B(3)*xdata_shifted); % Form of the equation
    % IC = [min(ydata) max(ydata) 0]; % Initial Conditions
    IC = [ic_baseline ic_c ic_m]; % Initial Conditions
    B = lsqcurvefit(F,IC,xdata_shifted,ydata,[lb_baseline lb_intersect lb_slope],[ub_baseline ub_intersect ub_slope], opts);

    Fitted_Contact_Point = (B(1) - B(2)) / B(3) + xdata_one; % in nm
    Baseline    = B(1);
    Intersect   = B(2);
    Slope       = B(3);  

    % Find closest point on curve to calculated contact point
    Crop_Length         = length(xdata_shifted);
    Shifted_Force_Curve = abs(Force_Curve(:,1) - Fitted_Contact_Point);
    Contact_Index       = find(Shifted_Force_Curve == min(Shifted_Force_Curve));
    Contact_Index       = Contact_Index(1);

    % Calculate where the contact point was found as a fraction of the length
    % of the cropped force curve. If the contact point was found in the
    % non-contact half of the force curve, the while loop will run again
    % (because Relative Position must be < 0.5, i.e. in the contact half).
    % Crop_Diff              = Data_Length - Crop_Length;
    % Relative_Contact_Point = Contact_Index - Crop_Diff;
    % Relative_Position_Crop = Relative_Contact_Point/Crop_Length;

    % Secondly, find out where the contact point was found as a fraction of the
    % length of the whole force curve. If it was found in the non-contact half,
    % the while loop goes again. This must be done in nm and not as an index
    % value, as in QNM the velocity of the tip is always changing.

    Data_Length_nm    = Force_Curve(end,1) - Force_Curve(1,1);
    Relative_Position = Force_Curve(Contact_Index,1)/Data_Length_nm;

    % If the contact point was not determined, more of the contact region is
    % included in the fitting of the Piecewise function, and the least squares
    % regression is attempted again.

    Attempt = Attempt + 1;

    Contact_Points(1,1) = Force_Curve(Contact_Index, 1);
    Contact_Points(1,2) = Force_Curve(Contact_Index, 2);


    Force_Curve_Corrected(:,1) = xTrace_Unfiltered - xTrace_Unfiltered(Contact_Index);
    Force_Curve_Corrected(:,2) = yTrace_Unfiltered - Baseline;

    % % Once the contact point and baseline have been determined, and the force
    % % curve has been corrected, the Z at minimum force is caluclated, which is
    % % the Z at contact. This is simply the Z difference, in nm, between Z at
    % % maximun force (which is given as the pixel value in the original height
    % % image) and Z at contact. These 'Zmin' values are later added to the
    % % 'Zmax' values to produce the 'True Height' image of the NPC (not within this function). 
    % Max_Indent_z_Index = find(Force_Curve_Corrected(:,2) == max(Force_Curve_Corrected(:,2))); %#ok<MXFND>
    % Max_Indent_z_Index = max(Max_Indent_z_Index); % if 2 of same force, takes the index of the point closest to the contact point
    % Indentation               = Force_Curve_Corrected(Max_Indent_z_Index,1);

    % Find the tip indentation in nm.
    Contact_Region_Z     = Force_Curve_Corrected(1:Contact_Index,1);
    Contact_Region_Force = Force_Curve_Corrected(1:Contact_Index,2);
    Max_Indent_z_Index   = find(Contact_Region_Force == max(Contact_Region_Force),1);
    Indentation          = abs(Contact_Region_Z(Max_Indent_z_Index));

    Contact_Points_Corrected(1,1) = Force_Curve_Corrected(Contact_Index,1);
    Contact_Points_Corrected(1,2) = Force_Curve_Corrected(Contact_Index,2);

end % end of while loop

end