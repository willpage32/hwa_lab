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

% Read in data from summary matricesU_plus_hre
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

fitfns = figure ; figure_format() ;
hold on ; plot(u_pr,V_pr) ; plot(u_po,V_po) ;
legend('Pre-experiment calibration','Post-experiment calibration')
title('Fitted calibration fns')

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

% surfer = figure;  figure_format(surfer) ;surf(T,V,u_fn)
% xlabel('time axis') ; ylabel('voltage axis') ; zlabel('velocity axis');
% title('Time linear interpolated U vs E correltation surface')

%% Load in the high reynolds data
[u_bar_hre,U_inf_hre,nu_hre,uvar_hre,x_hre,z_hre] = read_highRe();

% utau should be around 30 - n.hutch 445 
% upper limit at 15% of the layer thickness 
%% Read in the summary data
path_data = 'Data/';
[u_exp,v_exp,T_exp,T2_exp,P_exp,z_exp,f_exp] = read_summary_d(path_data);

%% Experimental constants 
k       = 0.4               ; % Von karman constant
A       = 5.0               ; % Clauser plot constant
R       = 286.9             ; % Individual Gas Const. -> ( J/ (kg K) )
T       = 20 + 273.15       ; % Temp in degrees Kelvin ()  
rho     = P_exp/(R*T)        % Fluid density -> (kg/m^3)
C1      = 1.458e-6          ; % Experimental? Constant from handout 
nu_air  = C1*(T^(1.5))/(rho*(T+110.5))  % 

%% Experimental clauser inputs 
close all

z_lre = (z_exp./1000) ; % Distance from wall (m) (assumed measured in mm) 
U_inf_lre = mean(u_exp) ;   %u_inf 
u_bar_lre = U_fnof_tnE((1:tf).',v_exp) ; % Corrleate voltage anemometor reading with velocity

y_plus_exp = (z_lre.*u_bar_lre)./(nu_air)   ; % EXPERIMENT Clauser x input
% scale by u_tau not u infinity, u+ should come to around ...
u_plusA_exp = u_bar_lre./U_inf_lre          ; % EXPERIMENT Clauser y input

figure ; semilogx(y_plus_exp , u_plusA_exp,'lineWidth',2) ;  hold on
% HIGH re
semilogx ((z_hre.*u_bar_hre)/nu_hre ,  u_bar_hre./U_inf_hre)

% figure_format(); axis([1.5e2,3e4,.4,1.01]) ;
title('Clauser Chart - Experimental data only ') ;
xlabel('y^{+} = (y u / \nu)') ; ylabel('U+ = u/u_{\infty}')

%% Guess U_tau

A=5;
k=k;

Cf = 0.000:.001:0.01

zplus_lre  = log(z_lre*U_inf_lre/nu_air);
zplus_hre  = log(z_hre*U_inf_hre/nu_air);

clauser_fn = @(Cf,zplus)   (1/k).*sqrt(Cf/2).*zplus + ...
              (1/k).*sqrt(Cf/2).*log(sqrt(Cf/2)) + ...
                   A*sqrt(Cf/2);
figure ; hold on ;

% size(u_hotwire/u_inf_lre)
% size(zplus_lre)
% size(clauser_fn(cf,zplus_lre))
plot(zplus_lre,u_bar_lre/U_inf_lre);
plot(zplus_lre,clauser_fn(Cf,zplus_lre),'--','color',[.7,.7,.7])
xlabel('z^+ for Low Reynold`s Number')
ylabel('u/U_{\infty}')
title('Clauser Chart for Low Reynold`s Number')
figure_format(); 


figure; hold on
plot(zplus_hre,u_bar_hre/U_inf_hre);
plot(zplus_hre,clauser_fn(Cf,zplus_hre),'--','color',[.7,.7,.7])
xlabel('z^+ for High Reynold`s Number')
ylabel('u/U_{\infty}')
title('Clauser Chart for High Reynold`s Number')
figure_format(); 


%%Old clauser stuff
% 
% clauser_x2 = exp( (1/k)  .* sqrt(Cf/2).*log(y.*u_inf./nu_air)  + ...
%         (1/k).*sqrt(Cf/2).*log(sqrt(Cf/2)) + A .* sqrt(Cf/2) );
% logx_clau2 = log(clauser_x2);
% hold on ;
% semilogx(clauser_x2,u_plusA_exp)
% break

% %% Find the HIGH RE friction velocity
% % High Reynolds data : [U_hre,Uinf_hre,nu_hre,uvar_hre,x_hre,z_hre] 
% 
% U_infplus_hre = U_hre/Uinf_hre      ; % First normalise by U_inftyf
% % cf_guess_hre = linspace(0.2,.5,5e2) ; % Vector of guesses for Cf
% cf_guess_hre = 0.4091;
% 
% y_plus_theory_hre = exp ((1/k) .* sqrt(cf_guess_hre/2).*log(z_hre.*Uinf_hre/nu_hre)  + ...
% 	(1/k).*sqrt(cf_guess_hre/2).*log(sqrt(cf_guess_hre/2)) + A * sqrt(cf_guess_hre/2) );
% 
% y_plus_hre = (z_hre.*U_hre)/nu_hre;
% 
% figure ; semilogx(y_plus_hre,U_infplus_hre) ; title('Experimental U+ , y+');
% xlabel('y+'); ylabel('U+');
% 
% figure ; semilogx(y_plus_theory_hre,U_infplus_hre,'-o') ; hold on ;
% semilogx(y_plus_hre,U_infplus_hre);
% title('Comparition of U+/y+ for varying cf to find U_{tau}');
% xlabel('y+'); ylabel('U+');
% 
% % 0.55 and 0.7 seem good U+ to seartch for linear regions
% 
% [mr_yplus_s, idx_yplus_st ] = min(abs(U_infplus_hre-0.55))
% [mr_yplus_e, idx_yplus_end] = min(abs(U_infplus_hre-0.7 ))
% 
% % Plot gradient evaluation points of experimental data on graph
% plot(y_plus_hre(idx_yplus_st) ,U_infplus_hre(idx_yplus_st) ,'k*','markerSize',10)
% plot(y_plus_hre(idx_yplus_end),U_infplus_hre(idx_yplus_end),'k*','markerSize',10)
% 
% % Gradient of experimental log region of velocity profile
% grad_hre_exp = (U_infplus_hre(idx_yplus_st) - U_infplus_hre(idx_yplus_end)) / ...
%     (log(y_plus_hre(idx_yplus_st)) - log(y_plus_hre(idx_yplus_end)))
% 
% % Find indicies of the start and end points of the linear region
% [mr_yplus_th_s, idx_yplus_th_st ] = min(abs(U_infplus_hre-0.55));
% [mr_yplus_th_e, idx_yplus_th_end] = min(abs(U_infplus_hre-0.7));
% 
% % Plot gradient evaluation points of all theory points on graph
% plot(y_plus_theory_hre(idx_yplus_st,:) ,U_infplus_hre(idx_yplus_st,:) ,'k*','markerSize',10)
% plot(y_plus_theory_hre(idx_yplus_end,:),U_infplus_hre(idx_yplus_end,:),'k*','markerSize',10)
% 
% % Evaluate all theoretical velcoity gradients 
% theory_hre_gradients = (U_infplus_hre(idx_yplus_th_end) - U_infplus_hre(idx_yplus_th_st) ) ./ ...
%     ( log(y_plus_theory_hre(idx_yplus_th_end,:)) - log(y_plus_theory_hre(idx_yplus_th_st,:)));
% 
% grad_diff_hre = grad_hre_exp-theory_hre_gradients;
% [mr_grad, idx_best_grad] = min(abs(grad_diff_hre));
% cf_found = cf_guess_hre(idx_best_grad);
% 
% % Plot the residual function for gradient difference
% figure ; plot(cf_guess_hre,-grad_hre_exp+theory_hre_gradients) ; hold on;
% title('Gradient difference - Residual plot');
% xlabel('Skin Friction - C_{f}') ; ylabel('Grandient difference function')
% y_dash_range = linspace(0.2,cf_found,2e2) ; y_dash_rangeB = zeros(size(y_dash_range));
% x_dash_range = linspace(-.01,0,2e2) ; x_dash_rangeB = cf_found*ones(size(x_dash_range));
% plot(y_dash_range,y_dash_rangeB,'--','color',[.7,.7,.7]) ;
% plot(x_dash_rangeB,x_dash_range,'--','color',[.7,.7,.7]) ;
% plot(cf_guess_hre(idx_best_grad),grad_diff_hre(idx_best_grad),'-p',...
%     'MarkerFaceColor',[255 105 180]./256,'MarkerSize',35);
% 
% % Create a new function to describe the matched cf y+ function
% y_plus_hre_theoryadj = exp( (1/k) .* sqrt(cf_found/2).*log(z_hre.*Uinf_hre/nu_hre)  + ...
% 	(1/k).*sqrt(cf_found/2).*log(sqrt(cf_found/2)) + A * sqrt(cf_found/2) );
% 
% % Plot the final guess for cf with the experimental data
% figure ; semilogx(y_plus_hre_theoryadj,U_infplus_hre) ; hold on ;
% semilogx(y_plus_hre,U_infplus_hre) ;
% title('Re-plot the experimental data with clauser theory, check gradients visually');
% 
% % Find Utau and re-normalise it
% rho_hre   = rho ; % Make the assumption both experiments are in air
% tau_w_hre = cf_found .*(1/2).*rho_hre.*(Uinf_hre.^2) ;  %[Pa]  or  N/m^2 
% U_tau_hre = sqrt(tau_w_hre./rho_hre) ;                  %[Pa]/[kg/m^3] = [m/s]
% Uplus_hre = U_hre/U_tau_hre ; % NORMALISE BY Utau
% 
% % Re-plot the velocity function in terms of u/utau
% figure ; semilogx(y_plus_hre, Uplus_hre)
% 
% %% find the gradients of clauser and experimental, find Cf
% 
% data_p_u = 21 ;
% data_p_l = 16 ;
% 
% x_th_1=logx_clau2(data_p_l,:);
% x_th_2=logx_clau2(data_p_u,:);
% 
% y_th_1=u_plusA_exp(data_p_l);
% y_th_2=u_plusA_exp(data_p_u);
% 
% x_ex_1 = logx_clau(data_p_l);
% x_ex_2 = logx_clau(data_p_u);
% 
% rise=y_th_2-y_th_1;
% run_th=x_th_2-x_th_1;
% run_ex=x_ex_2-x_ex_1;
% 
% grad_th=rise./run_th;
% grad_ex=rise./run_ex;
% 
% grad_diff=grad_th-grad_ex;
% 
% [e Cf_loc]=min(abs(grad_diff));
% Cf_true=Cf(Cf_loc)
% 
% % figure ; semilogx()
% 
% figure;
% hold on
% plot(Cf,grad_diff,'*')
% plot(Cf,zeros(size(Cf)));
% plot(Cf(Cf_loc),grad_diff(Cf_loc),'-p',...
%     'MarkerFaceColor',[255 105 180]./256,...
%     'MarkerSize',35)
% title('Difference Between Log Region Gradients for Varying C_f');
% xlabel('C_f')
% ylabel('Difference (unitless)')
% % figure_format(); 
% 
% %% Caluclate Tau and U_tau
% % once Cf is found out, we can calculate tau_w and then U_tau
% tau_w = Cf_true .*(1/2).*rho.*(u_inf_lre.^2) ;  % [Pa]  or  N/m^2 
% U_tau = mean(sqrt(tau_w./rho)) ;                  % [Pa]/[kg/m^3] = [m/s]
% U_exp_u_tau = u_hotwire/U_tau ;             % 
% 
% figure ; semilogx(y_plus_exp,U_exp_u_tau) ;
% %% Define true Clauser plot
% 
% clauser_true =(1/k) * sqrt(Cf_true/2).*log(z_lre.*u_inf_lre/nu_air)  + ...
%         (1/k).*sqrt(Cf_true/2).*log(sqrt(Cf_true/2)) + A * sqrt(Cf_true/2);
% 
% %% Variance
% 
% %load('times.mat')
% load('volts.mat')
% 
% voltage_variance=var(volt_full);
% U_var = U_fnof_tnE((1:tf).',voltage_variance') ;
% 
% %%  Plot Variance QN3/4
% figure; hold on
% plot(z_lre,U_var)  % what units are zexp and zhre in?
% plot(z_hre,uvar_hre)
% xlabel('Distance From Wall [m]')
% ylabel('Variance')
% 
% figure;
% subplot(2,1,1)
% plot(z_lre,U_var)  % what units are zexp and zhre in?
% ylabel('Variance')
% subplot(2,1,2)
% plot(z_hre,uvar_hre)
% ylabel('Variance')
% xlabel('Distance From Wall [m]')
% 
% % 1/20 for low re, 1/40 for high re.
% % 2 methods are eqiuvlant.
% 
% %% Plot Clauser Graphs (QN 1 & 2)
% 
% % % HIGH re
% % figure ; semilogx ((z_hre.*U_hre)/nu_hre ,  U_hre./Uinf_hre)
% % hold on
% % title('Clauser plot, high re data') ;
% % legend('High Re')
% 
% figure ;    
% semilogx(exp(clauser_true),u_plusA_exp) ;
% title('THEORETICAL Inner scaled velocity variance') ; xlabel('y+') ; ylabel('U+') ;
% axis([5e3,1e6,0.4,1.02]) ;
% 
% figure ; 
% deficit_y = (u_inf_lre-u_hotwire)./(U_tau) ;
% semilogx(exp(clauser_true),deficit_y)  ;
% title('THEORETICAL Outer scaled velocity profile (deficit form)');
% ylabel('u-uinf/utau') ; xlabel('y+') ;
% 
% %semilogx()
% % shift = 3.3 ;
% % semilogx(exp(log(clauser_x)+shift),clauser_y) ;
% % legend('Clauser Theory','Experimental') ;
% % xlabel('y+') ; ylabel('U+') ; 
% 
% % % Find experimental constant 
% % Uplus_known = u_hotwire(data_p_l) / u_inf(data_p_l)
% % yplus_known = (y(data_p_l)*u_hotwire(data_p_l)) / ( 1e2*nu_air)
% % 
% % log_range = data_p_l:data_p_u
% % logbit_uplus = (u_hotwire(log_range)) ./ (u_inf(log_range));
% % logbit_yplus = (y(log_range).*u_hotwire(log_range))/ nu_air;
% % 
% % [xData, yData] = prepareCurveData( logbit_yplus, logbit_uplus );
% % ft = fittype( 'a + .1384*log(x)', 'independent', 'x', 'dependent', 'y' );
% % [log_bit_fit, gof] = fit( xData, yData, ft );
% % 
% % rng = linspace(.1,3e4,1e2);
% % figure ; plot(clauser_x,clauser_y) ; hold on ;
% % plot(rng,log_bit_fit(rng))
% % 
% % figure ; semilogx(rng,log_bit_fit(rng)) ; hold on ;
% % semilogx(clauser_x,clauser_y)
% 
% % plot( log_bit_fit, xData, yData );
% % semilogx(log_bit_fit())
% % syms const1
% % 
% % ayy = solve(Uplus_known==grad_ex*yplus_known + const1)
% % ayy = vpa(ayy,3)
% % 
% % fplot(@(yplus) .1384*log(yplus)+ayy ,[1e1,1e6])
% 
% % x2u = u_hotwire(data_p_u)%=21;
% % x1u = u_hotwire(data_p_l)%=16;
% % 
% % x1 = log((y(data_p_l)*x1u)/nu_air); 
% % y1 = x1u ./ u_inf(data_p_l) ;
% % x2 = log((y(data_p_u)*x2u)/nu_air); 
% % y2 = x2u ./ u_inf(data_p_u);
% % 
% % [p,~,mu] = polyfit([x1, x2], [y1, y2], 1);
% % % a = coefficients(1);
% % % b = coefficients(2);
% % f = polyval(p,[0:.001:1],[],mu);
% % plot([0:.001:1],f)
% % 
% % [acor,lag] =xcorr(clauser_true,log(clauser_x));
% % [acor,lag']
% 
% %% Qn6.7a, determine \delta_{99}, \delta*, \theta , H , Cf
% %Determine for low Re (experimental data)
% 
% %Delta 99 - Bpoundary Layer thickness
% vel_ratio=u_hotwire./u_inf_lre;
% [e loc99]=min(abs(vel_ratio-.99));
% delta99=z_lre(loc99)
% 
% %show the velocity on a plot
% figure; hold on; xlabel('Distance from Wall [m]');ylabel('Velocity')
% plot(z_lre,u_hotwire); plot(delta99,u_hotwire(loc99),'r*')
% legend('Hotwire Velocity','99% U_0'); 
% 
% %Delta* - Displacement Thickness
% %uses trapz to integrate w.r.t y
% T1=(1-u_hotwire./u_inf_lre);
% delta_star=trapz(z_lre,T1);
% 
% %\Theta - Momentum thickness
% T2=(u_hotwire./u_inf_lre).*(1-u_hotwire./u_inf_lre);
% Theta=trapz(z_lre,T2);
% 
% %H - Shape factor
% %The higher the value of H, the stronger the adverse pressure gradient. 
% %A high adverse pressure gradient can greatly reduce the Reynolds number at 
% %which transition into turbulence may occur.
% H=delta_star/Theta;
% 
% %Skin Friction Coefficient
% Cf_true;
% 
% %% Qn6.7b, determine for high Re data 
% %[U_hre,Uinf_hre,nu_hre,uvar_hre,x_hre,z_hre] = read_highRe();
% nu_air./U_tau
% 
% %Delta 99 - Bpoundary Layer thickness
% vel_ratio_hre=U_hre./Uinf_hre;
% [e_hre loc99_hre]=min(abs(vel_ratio_hre-.99));
% delta99_hre=z_hre(loc99_hre)
% 
% %show the velocity on a plot
% figure; hold on; xlabel('Distance from Wall [m]');ylabel('Velocity')
% plot(z_hre,U_hre); plot(delta99_hre,U_hre(loc99_hre),'r*')
% legend('Hotwire Velocity High Re','99% U_0'); 
% 
% %Delta* - Displacement Thickness
% %uses trapz to integrate w.r.t y
% T1_hre=(1-U_hre./Uinf_hre);
% delta_star_hre=trapz(z_hre,T1_hre);
% 
% %\Theta - Momentum thickness
% T2_hre=(U_hre./Uinf_hre).*(1-U_hre./Uinf_hre);
% Theta_hre=trapz(z_hre,T2_hre);
% 
% %H - Shape factor
% %The higher the value of H, the stronger the adverse pressure gradient. 
% %A high adverse pressure gradient can greatly reduce the Reynolds number at 
% %which transition into turbulence may occur.
% H=delta_star_hre/Theta_hre;
% 
% %Skin Friction Coefficient
% Cf_true;
% 
% retau= delta99_hre*U_tau/nu_air
% 
% %% Qn 7
% 
% %find the data point
% [e loc90]=min(abs(vel_ratio-.90));
% 
% %read in the daq
% data_string=['Data/'];
% daq_string=['.daq'];
% data_n=num2str(loc90,'%i');
% data_loc = strcat(data_string,data_n,daq_string);
% 
% 
% [daq_90,time,abstime] = daqread(data_loc);
% 
% mean_90=mean(daq_90);
% mean_90=mean(daq_90);
% histogram(daq_90(:,1))


% %% Functions
% 
% function output = clauser_fn(Cf)
%               (1/k).*sqrt(Cf/2).*z_plus + ...
%               (1/k).*sqrt(Cf/2).*log(sqrt(Cf/2)) + ...
%                    A*sqrt(Cf/2);end