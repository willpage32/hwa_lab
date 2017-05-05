% William Page (587000) - Reads the summary matrix from specified path 

function [u1,v1,T1,T2,P1,z1,f1] = read_summary_d(path)

file_call = [path,'/summary.mat'] ; % Define path for data
load(file_call)                   ; % Load it in 

% Assign variables
u1 = u; v1 = volt; T1 = T_i ; T2 = T_o ; P1 = P_atm ; f1 = freq ; z1 = z_meas ;

% Clip weird zeros
i_zeros_u = find(u1 == 0)           ; % Find the zeros (except the first datapoint)
i_zeros_u = i_zeros_u(i_zeros_u~=1) ; % Take our any data in element 1
u1(i_zeros_u) = []                  ; % Replace them with empty

i_zeros_v = find(v1 == 0)           ; % Find where volts are == 0  
i_zeros_v = i_zeros_v(i_zeros_v~=1) ; % Take our any data in element 1
v1(i_zeros_v) = []                  ; % Replace with empty