% William Page (587000) - Kevin Rassool (540733)   ;
% Semester 2 2017 - University of Melbourne        ; Started:     01/5/17
% MCEN90018 - Advanced Fluid Dynamics              ; Last Edited: 19/5/17
% Hot Wire Laboratory : Turbulent behaviours
%
% Main Script

%% Read in data
clc, clear , close all

% Define data directories
path_pre  = 'Data/pre'  ; % Pre-calibration path 
path_post = 'Data/post' ; % Post-calibration path 

% Read in data from summary matrices
[u_pre ,V_pre ,T_pre ,P ]         = read_summary(path_pre)  ; % Pre -calibration  
[u_post,V_post,T_post,P_atm_post] = read_summary(path_post) ; % Post-calibration

%% Define pre and post fits
% Make a fit, fitType == 1 for kings law , == 2 for 3rd order poly
close all
fitType = 2     ; % 1 = Kings Law ; 2 = 3rd Order poly fit

%  Make a kings law / polyn fit 
V_pr_fit = HWA_Calib_fitter(u_pre ,V_pre , fitType) ;
V_po_fit = HWA_Calib_fitter(u_post,V_post, fitType) ;

u_pr_fitdat = linspace(min(u_pre) , max(u_pre) , 1e2) ;
u_po_fitdat = linspace(min(u_post), max(u_post), 1e2) ;

V_pr_fitdat = V_pr_fit(u_pr_fitdat) ;
V_po_fitdat = V_po_fit(u_po_fitdat) ;

% Plot of pre calibration, OG data with curve through points
figure ; hold on ; plot(V_pr_fitdat,u_pr_fitdat)  ;
plot(V_pre,u_pre,'o') ; title('Calibration plot, pre experiment data')
xlabel('Voltage (V)') ; ylabel('Flow Speed (m/s)') ;
legend('data-fit','original data') ;

% Plot of post calibration, OG data with curve through points
figure ; hold on ; plot(V_pr_fitdat,u_po_fitdat) ;
plot(V_post,u_post,'o') ; title('Calibration plot, post experiment data')
xlabel('Voltage (V)') ; ylabel('Flow Speed (m/s)') ;
legend('data-fit','original data') ;

% Plot the pre and post calibration signals on the same axes
figure ; hold on ; plot(V_pr_fitdat,u_pr_fitdat)  ;
plot(V_po_fitdat,u_po_fitdat) ; title('Fitted calibration curves')
xlabel('Voltage (V)') ; ylabel('Flow Speed (m/s)') ;
legend('pre-experiment','post-experiment') ;

% for a fit at time 't' during the experiment
t0 = 0  ; % Start time (mins)
tf = 30 ; % End time   (mins)
time_start = t0*ones(length(V_pr_fitdat),1);
time_end   = tf*ones(length(V_po_fitdat),1);

V_fit = [V_pr_fitdat;V_po_fitdat]    ; % Voltage  fit functions
u_fit = [u_pr_fitdat.';u_po_fitdat.']; % Velocity fit functions
t_fit = [time_start;time_end]        ; % Time vector to match input curves

%% Make a polyn fit 
[xd1,yd1] = prepareCurveData(V_pre ,u_pre ) ; % Prepare data, strips inf/nan/etc
[xd2,yd2] = prepareCurveData(V_post,u_post) ; % Prepare data, strips inf/nan/etc
[V_prfit, ~] = fit(xd1, yd1, fittype('poly3')) ; % Create the  pre-fit
[V_pofit, ~] = fit(xd2, yd2, fittype('poly3')) ; % Create the post-fit

% Velocity fucntions describing some new 'u'
u_pr = linspace(min(u_pre) ,max(u_pre) ,1e2).' ; 
u_po = linspace(min(u_post),max(u_post),1e2).' ; 
V_pr = V_prfit(u_pr) ; V_po = V_pofit(u_po)  ;  

figure ; hold on ; plot(u_pr,V_pr) ; plot(u_po,V_po) ; title('Fitted calibration fns')
onesies = ones(length(u_pr),1);
t0 = 0 ; tf = 30 ; 
t0s = t0*onesies ; tfs = tf*onesies ;

% Define column vectors to create my surface function for U vs E
u_fit = [u_pr;u_po];
E_fit = [V_pr;V_po];
t_fit = [t0s ; tfs];

[xData, yData, zData] = prepareSurfaceData( E_fit, u_fit, t_fit );
ft = fittype( 'poly22' );% Set up fittype and options.
[fitresult, gof] = fit( [xData, yData], zData, ft );

%% Clauser plot input

% Intialise matrix of daqs
ndaqs   = 30      ;
ndaqpts = 30*30e3 ;
data_full = zeros(ndaqs*ndaqpts,2);
for i=1:ndaqs
    filepath  = ['Data/',num2str(i),'.daq'];
    data_temp = daqread(filepath);
    strow     = (i-1)*(ndaqpts)+1;
    endrow    = strow + ndaqpts-1;
    data_full(strow:endrow,:) = data_temp ;
end

figure(plot())

mu = 1.81e5 ;
nu = mu/rho ;

% U    = ?
% Uinf = ?
% log(y*U_inf/mu)

figure ; plot(data_full(:,2),data_full(:,1))

%% Load in the high reynolds data
[U_hre,Uinf_hre,nu_hre,uvar_hre,x_hre,z_hre] = read_highRe();

R   = 8.314459848 ; % Gas Const. -> ( kg m^(2) / ( s^2 K mol) )
T   = 20 + 273.15 ; % Temp in degrees celcius -> (degC) 
rho = P/(R*T)     ; % Fluid density -> (kg/m^3)

figure ; hold on ;
plot(log((z_hre.*U_hre)/nu_hre) ,  U_hre./Uinf_hre)