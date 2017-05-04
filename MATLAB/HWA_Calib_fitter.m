% William Page (587000) - Kevin Rassool (540733)   ;
% Semester 2 2017 - University of Melbourne        ; Started:     01/5/17
% MCEN90018 - Advanced Fluid Dynamics              ; Last Edited: 19/5/17
% Hot Wire Laboratory : Turbulent behaviours
%
% Curve Fitting: Used for pre-cal and post cal
%
% FitType == 1 Provides a kings law fit
% FitType == 2 Provides a 3rd order poly fit
%
% Inputs are : HWA_Calib_fitter(VELOCITY,VOLTAGE,fitType) 

function V_handle = HWA_Calib_fitter(u,V,fitType) 

% Make a custom fit
[xd1,yd1] = prepareCurveData(u,V) ; % Prepare data, strips inf/nan/etc
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );

% Specify fit type
if fitType==1
    ft = fittype( '( app + bpp*x^(1/n) )', 'independent', 'x', 'dependent', 'y' );
    opts.startPoint = [] % First guess for fit co-efficients
elseif fitType == 2
    ft = fittype( 'a*x^3 + b*x^2 c*x + d', 'independent', 'x', 'dependent', 'y' );
    opts.startPoint = [] % First guess for fit co-efficients
else
    fprintf('Wrong input')
end



% Make fits 
[V_handle, ~] = fit( xd1, yd1 , ft ,opts);
end