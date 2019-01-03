% This function applies the Hertz model to a force curve. The force curve
% should be already corrected for the baseline and contact point, and
% should be in the format nN vs nm (tip-sample separation). The input
% arguments for the function are: the corrected force curve, entered as an
% Nx2 array; the index for its contact point along the x axis (tip-sample
% separation); and the length of indentation, in nm, that is to be
% modelled. Typically, 20 nm is modelled. The output is Young's Modulus
% in MPa.

function [YM] = HertzModel(Force_Curve_Corrected, Contact_Index, Indentation_For_Fit)

%%%======Choose indentation for fitting Hertz model(nm)======%%%
Indentation_Hertz_Fit = Indentation_For_Fit;
%%%==========================================================%%%

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