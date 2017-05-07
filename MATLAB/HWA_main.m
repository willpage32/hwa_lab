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
[u_pre ,V_pre ,Tpr	 ,Ppr ]       = read_summary(path_pre)  ; % Pre -calibration  
[u_post,V_post,T_post,P_atm_post] = read_summary(path_post) ; % Post-calibration

%% Make a polyn fit 
[xd1,yd1] = prepareCurveData(V_pre , u_pre ) ; % Prepare data, strips inf/nan/etc
[xd2,yd2] = prepareCurveData(V_post, u_post) ; % Prepare data, strips inf/nan/etc
[u_prfit, Gpr] = fit(xd1, yd1, fittype('poly3')) ; % Create the  pre-fit
[u_pofit, Gpo] = fit(xd2, yd2, fittype('poly3')) ; % Create the post-fit

% Velocity fucntions describing some new 'u'
V_pr = linspace(min(u_pre) ,max(u_pre) ,1e2).' ; 
V_po = linspace(min(u_post),max(u_post),1e2).' ; 

u_pr = u_prfit(V_pr) ; u_po = u_pofit(V_po)  ;  

figure ; hold on ; plot(u_pr,V_pr) ; plot(u_po,V_po) ; title('Fitted calibration fns')
legend('Pre-experiment calibration','Post-experiment calibration')
onesies = ones(length(V_pr),1);
t0 = 0 ; tf = 30 ; 
t0s = t0*onesies ; tfs = tf*onesies ;

%% Fit a surface to the data so we can evaluate at any time
u_fit = [u_pr;u_po] ; E_fit = [V_pr;V_po] ; t_fit = [t0s ; tfs] ;
[xData, yData, zData] = prepareSurfaceData( t_fit, E_fit, u_fit );
[U_fnof_tnE, ~] = fit( [xData, yData], zData, fittype('poly15'));

time_fn = linspace(t0,tf,1e2).';
volt_fn = linspace(min(u_fit),max(u_fit),1e2).';

[T,V] = meshgrid(time_fn,volt_fn);
u_fn    = U_fnof_tnE(T,V);

figure ; plot(U_fnof_tnE,[t_fit,E_fit],u_fit) ;
xlabel('time axis') ; ylabel('voltage axis') ; zlabel('velocity axis');
title('Time linear interpolated U vs E correltation surface inputs')

% figure; surf(T,V,u_fn)
% xlabel('time axis') ; ylabel('voltage axis') ; zlabel('velocity axis');
% title('Time linear interpolated U vs E correltation surface')

%% Load in the high reynolds data
[U_hre,Uinf_hre,nu_hre,uvar_hre,x_hre,z_hre] = read_highRe();

figure ; plot(log((z_hre.*U_hre)/nu_hre) ,  U_hre./Uinf_hre)
title('Clauser plot, high re data')

% Read in the esummary data
path_data = 'Data/';
[u_exp,v_exp,T_exp,T2_exp,P_exp,z_exp,f_exp] = read_summary_d(path_data);

% Experimental constants 
R       = 8.314459848       ; % Gas Const. -> ( kg m^(2) / ( s^2 K mol) )
T       = 20 + 273.15       ; % Temp in degrees celcius -> (degC) 
rho_hre = P_exp/(R*T_exp)   ; % Fluid density -> (kg/m^3)
nu_air  = 15.11e-6          ; % Dynamic viscoity ( m^2 / s) 

%% Experimental clauser inputs 
clasuer_x = log( ((z_exp./1000).*u_exp)./(nu_air) ) ; % Clauser x input
% clauser_y = u_exp/u_exp_inf;

figure; plot(u_exp) ; title('Averaged u velocity') ; 
xlabel('Currently n (-) should be distance (m)') ; ylabel('U/U_{\inf} (-)')

%% Use of the clauser plot 
