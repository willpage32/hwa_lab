%% Calibration - Exponentials
close all
% Fit exponential functions for the data
V_u_pre_fit  = fit(V_pre, u_pre, 'exp1') ; % Exponential fit for pre data
V_u_post_fit = fit(V_post,u_post,'exp1') ; % Exponential fit for post data

% Fit piecewise linear functions to compare signals
fit_vpre_lin = fit(V_pre, u_pre,'linearinterp') ; 

% Define vectors to evalute fits at 
V_prefit  = linspace(min(V_pre) ,max(V_pre) ,1e2); % Fit u data to eval Vpre  at
V_postfit = linspace(min(V_post),max(V_post),1e2); % Fit u data to eval Vpost at

% Plot original functions
figure ; hold on    ;
plot(V_pre,u_pre)   ; 
plot(V_post,u_post) ; axis equal tight

figure ;
hold on ;
plot(V_prefit,V_u_pre_fit(V_prefit)); 
plot(V_postfit,V_u_post_fit(V_postfit)) ; axis equal tight ; 

dVpre = V_u_pre_fit(V_prefit) - fit_vpre_lin(V_prefit);
% Plot fit functions
% Use polyfit, order 'n' poly with co-efs
% drp = linspace(0,13,1e2+1);
% u_volt_pre_coeff = polyfit(u_pre,volt_pre,3);
% fit_pre = polyval(u_volt_pre_coeff,drp);

% u_volt_pre_fitfn1 = fit(u_pre,volt_pre,'poly1');

% plot(u_volt_pre_fitfn2,u_pre,volt_pre)
% Interpolated function 