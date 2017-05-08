% William Page (587000) - Kevin Rassool (540733)   ;
% Semester 2 2017 - University of Melbourne        ; Started:     01/5/17
% MCEN90018 - Advanced Fluid Dynamics              ; Last Edited: 19/5/17
% Hot Wire Laboratory : Turbulent behaviours
%
% High Frequency Data Analysis Script

clc, clear, close all   % Clear the MATLAB envirnoment

ndaqs       = 30      ; % Number of data files (maybe 30 becuase 30 seconds?)
ndaqpts     = 30*30e3 ; % Number of points in each .daq file
data_full   = zeros(ndaqs*ndaqpts,2) ; % Intialise matrix of daqs

for i=1:ndaqs % Read all high frequency data into one vector
    filepath  = ['Data/',num2str(i),'.daq'] ; % File path for Daq file
    data_temp = daqread(filepath) ; % Temporary variable for data
    strow     = (i-1)*(ndaqpts)+1 ; % Starting row number for array 
    endrow    = strow + ndaqpts-1 ; % Ending row number for array
    data_full(strow:endrow,:) = data_temp ; % Complete data set
end
%%
close all ;

% Define experiment constants
R       = 8.314459848       ; % Gas Const. -> ( kg m^(2) / ( s^2 K mol) )
Temp    = 20 + 273.15       ; % Temp in degrees celcius -> (degC) 
P_air   = 101325            ; % Ambient air pressure 
rho = P_air/(R*Temp)        ; % Fluid density -> (kg/m^3)
mu  = 1.81e5 ; nu = mu/rho  ; % Viscocities

freq    = 30e3 ; % Data sampling frequency (for HWA!, not dyn pressure)
up2time = .01 ; % Time up to take the data 
up2n = freq*up2time  % Number of elements to pull out of the arrays
clip_1 = data_full(1:up2n,1) ; % First  row of clipped data
clip_2 = data_full(1:up2n,2) ; % Second row of clipped data

figure ; plot(clip_1) ; 
title('.daq file data, column 1, (HWA Voltage I think?)') ;
xlabel('Time (1/30000 s)') ; ylabel('Voltage (V)') ;

figure ; plot(clip_2) ;
title('.daq file data, column 2,first 200 points') ;
xlabel('Time (1/30000 s)') ; ylabel('Dynamic pressure') ;

%% Energy spectra analysis

G1 = fft(clip_1);
A1 = sqrt(4*(G1./up2n).*conj(G1./up2n));

leq_nyq = 2:1:((up2n/2)+1);

figure ; plot(A1(leq_nyq)) ; title('Amplitude of fourier modes, sub nyquist mode')
xlabel('Sample number') ; ylabel('Mode amplitude') ;

figure ; plot(A1(2:(end-1))) ; title('ALL fourier mode amplitudes except the mean')

%  xlabel('Sample Number') ; ylabel('who knows') ;
% legend('Row 1','Row 2')

% U    = ?
% Uinf = ?
% log(y*U_inf/mu)


