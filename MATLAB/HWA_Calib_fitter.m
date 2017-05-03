% William Page (587000) - Kevin Rassool (540733)   ;
% Semester 2 2017 - University of Melbourne        ; Started:     01/5/17
% MCEN90018 - Advanced Fluid Dynamics              ; Last Edited: 19/5/17
% Hot Wire Laboratory : Turbulent behaviours
% Curve Fitting: Using an exponential, 3rd order poly or sqrt fn
%
% FitType == 1 Provides a sqrt fit
% FitType == 2 Provides a 3rd order poly fit
%
% plotstho == 1 provides plots for the fits
function V_handle = HWA_Calib_fitter(u,V,fitType,plotstho) 

% Make a custom fit
[xd1,yd1] = prepareCurveData(u ,V ) ; % Prepare data, strips inf/nan/etc

% Specify fit type
if fitType==1
    ft = fittype( 'a + b*((x)^(1/2))', 'independent', 'x', 'dependent', 'y' );
elseif fitType == 2
    ft = fittype( 'a*x^3 + b*x^2 c*x + d', 'independent', 'x', 'dependent', 'y' );
else
    fprintf('Wrong input')
end

opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.startPoint = [u(1) V(1)]; % Make the first point extra correct

% Make fits 
[V_handle, ~] = fit( xd1, yd1 , ft, opts );

% Plot the fit if asked to
if strcmp(plotstho ,'yeahm8' )
    u_fitdat = linspace(min(u), max(u), 1e2);
    figure ; hold on ; 
    plot(u,V); axis equal;
    plot(u_fitdat,V_handle(u_fitdat)) ;
    axis([0,14,-6,8]); 
    title('Data from summary matrices') ; 
    xlabel('u velocity') ; ylabel('Voltage');
    legend('Original Data','Curve Fit');
elseif strcmp(plotstho , 'nah')
    
end