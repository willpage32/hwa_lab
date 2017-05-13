% William Page (587000) - Kevin Rassool (540733)   ;
% Semester 2 2017 - University of Melbourne        ; Started:     01/5/17
% MCEN90018 - Advanced Fluid Dynamics              ; Last Edited: 19/5/17
% Hot Wire Laboratory : Turbulent behaviours
%
% High Frequency Data Analysis Script

clc, clear, close all   % Clear the MATLAB envirnoment

ndaqs       = 30      ; % Number of data files (maybe 30 becuase 30 seconds?)
ndaqpts     = 30*30e3 ; % Number of points in each .daq file
% data_full   = zeros(ndaqs*ndaqpts,2) ; % Intialise matrix of daqs
[time_full,volt_full] = deal(zeros(ndaqpts,ndaqs));


for i=1:ndaqs % Read all high frequency data into one vector
    filepath  = ['Data/',num2str(i),'.daq'] ; % File path for Daq file
    [data_temp,time] = daqread(filepath) ; % Temporary variable for data
    time_full(:,i) = data_temp(:,1) ; % Build the time data
    volt_full(:,i) = data_temp(:,2) ; % Build the voltage data
end
%%
save('time','times_full')
save('volt','volts_full')
%%
close all ;

% Define experiment constants
R       = 8.314459848       ; % Gas Const. -> ( kg m^(2) / ( s^2 K mol) )
Temp    = 20 + 273.15       ; % Temp in degrees celcius -> (degC) 
P_air   = 101325            ; % Ambient air pressure 
rho = P_air/(R*Temp)        ; % Fluid density -> (kg/m^3)
mu  = 1.81e5 ; nu = mu/rho  ; % Viscocities

freq = 30e3 ; % Data sampling frequency (for HWA!, not dyn pressure)
dt   = .01  ; % Time up to take the data 
N    = freq*dt ; % Number of elements to pull out of the arrays
% t_st = 1       ; % Start time to take samples from
% t_f  = s_st+N*10;

dr   = 1:N ; % Data range to take form master array

clip_1 = data_full(dr,1) ; % First  row of clipped data
clip_2 = data_full(dr,2) ; % Second row of clipped data

figure ; plot(clip_1) ; 
title('.daq file data, column 1, (HWA Voltage I think?)') ;
xlabel('Time (1/30000 s)') ; ylabel('Voltage (V)') ;
axis([0,300,.51,.55]);

figure ; plot(clip_2) ;
title('.daq file data, column 2,first 200 points') ;
xlabel('Time (1/30000 s)') ; ylabel('Dynamic pressure') ;

figure ; 
for pv = 1:100
    st      = pv*N;
    endr    = st+N;
    clip_r  = st:endr;
    
    plot(data_full(clip_r))
    axis([0,300,.51,.55]);
    pause
end

%% Energy spectra analysis

G1 = fft(clip_1);
A1 = sqrt(4*(G1./N).*conj(G1./N));

leq_nyq = 2:1:((N/2)+1);

figure ; plot(A1(leq_nyq)) ; title('Amplitude of fourier modes, sub nyquist mode')
xlabel('Sample number') ; ylabel('Mode amplitude') ;

figure ; plot(A1(2:(end-1))) ; title('ALL fourier mode amplitudes except the mean')

%  xlabel('Sample Number') ; ylabel('who knows') ;
% legend('Row 1','Row 2')

% U    = ?
% Uinf = ?
% log(y*U_inf/mu)


