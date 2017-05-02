% William Page (587000) - Kevin Rassool (xxxxxx) ;
% Semester 2 2017 - University of Melbourne        ; Started:     01/5/17
% MCEN90018 - Advanced Fluid Dynamics              ; Last Edited: 19/5/17
% Hot Wire Laboratory : Turbulent behaviours
%
% Main Script

clc, clear, close all

load('highRE.mat') % Load high Re number data from LMS

whos

% highRE LMS variables
U       = U    ; % Velocity at some z       (m/s)
Uinf    = Uinf ; % Free stream velocity     (m/s)
nu      = nu   ; % Kinematic viscosity      (mu/rho) = (m^2/s)
uvar    = uvar ; % Variance of the signal U (-)
x       = x    ; % Single numer == 21 ?     (m?)
z       = z    ; % Z height                 (m)

