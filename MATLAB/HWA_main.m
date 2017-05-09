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
title('Clauser plot, high re data') ;

% Read in the esummary data
path_data = 'Data/';
[u_exp,v_exp,T_exp,T2_exp,P_exp,z_exp,f_exp] = read_summary_d(path_data);

% Experimental constants 
k       = 0.4               ; % Von karman constant
A       = 5.0               ; % Clauser plot constant
R       = 286.9             ; % Individual Gas Const. -> ( J/ (kg K) )
T       = 20 + 273.15       ; % Temp in degrees Kelvin ()  
rho     = P_exp/(R*T)       ; % Fluid density -> (kg/m^3)
nu_air  = 15.11e-6          ; % Dynamic viscoity ( m^2 / s) 
%% Experimental clauser inputs 
close all

y = (z_exp./1000); % Distance from wall (m) (assumed

u_hotwire = U_fnof_tnE((1:tf).',v_exp) ;
clauser_x = (y.*u_hotwire)./(nu_air)   ; % Clauser x input
u_inf     = u_exp ;
clauser_y = u_hotwire./u_inf ;

logxc = log(clauser_x);

figure ; plot(logxc,clauser_y) ; title('experimental Clauser plot') ; 
ylabel('U^{+}') ; xlabel('y^{+}') ; hold on ;

Cf = .3741  ;
clauser_x2 =(1/k) * sqrt(Cf/2)*log(y.*u_inf/nu_air)  + ...
        (1/k)*sqrt(Cf/2)*log(sqrt(Cf/2)) + A * sqrt(Cf/2) ;
plot(clauser_x2,clauser_y); 

% rise =
% run  = 
% grad_experi =rise 
% grad_theory = 

% set(gca,'xscale','log')
% once Cf is found out, we can calculate tau_w and then U_tau
tau_w = Cf*(1/2)*rho*(u_inf.^2) ;
U_tau = sqrt(tau_w/rho) ;
