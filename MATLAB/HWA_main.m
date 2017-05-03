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
% Make a fit, fitType == 1 for sqrt, ==2 for 3rd order poly

plotstho = 'nah'; % Plots the fitted curve vs orginal for yes, = 'yeahm8'
fitType = 2     ; % Make a sqrt fit, use =2 for 3rd order poly fit

V_pre_fit  = HWA_Calib_fitter(u_pre, V_pre , 1,plotstho);
V_post_fit = HWA_Calib_fitter(u_post,V_post, 1,plotstho);

u_pr_fitdat = linspace(min(u_pre) , max(u_pre) , 1e2);
u_po_fitdat = linspace(min(u_post), max(u_post), 1e2);

figure ; hold on ;
plot(u_pr_fitdat,V_pre_fit(u_pr_fitdat))  ;
plot(u_po_fitdat,V_post_fit(u_po_fitdat)) ; hold off ;

mu = 1.81e5;
nu = mu/rho;

figure ; plot()

%% Load in the high reynolds data
[U_hre,Uinf_hre,nu_hre,uvar_hre,x_hre,z_hre] = read_highRe();

R   = 8.314459848 ; % Gas Const. -> ( kg m^(2) / ( s^2 K mol) )
T   = 20 + 273.15 ; % Temp in degrees celcius -> (degC) 
rho = P/(R*T)     ; % Fluid density -> (kg/m^3)

figure ; hold on ;
plot(log((z_hre.*U_hre)/nu_hre) ,  U_hre./Uinf_hre)

