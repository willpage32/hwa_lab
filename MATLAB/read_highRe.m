% William Page (587000) - Reads the high re data into good varnames

function [U_hre,Uinf_hre,nu_hre,uvar_hre,x_hre,z_hre] ... 
            = read_highRe()
load('highRE.mat')
% HighRE LMS variables
U_hre    = U    ; % Velocity at some z       (m/s)
Uinf_hre = Uinf ; % Free stream velocity     (m/s)
nu_hre   = nu   ; % Kinematic viscosity      (mu/rho) = (m^2/s)
uvar_hre = uvar ; % Variance of the signal U (-)
x_hre    = x    ; % Single numer == 21 ?     (m?)
z_hre    = z    ; % Z height                 (m)
