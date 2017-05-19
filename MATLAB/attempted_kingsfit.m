% Define pre and post fits
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