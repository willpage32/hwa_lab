




%u_bar=data_in
%U1=data_in     %U_inf

%universal constants
k=0.4;
A=5;
%properties of air
ro=1.205        %kg/m^3 density
mu=15.11e-6     %m^2/s  kinematic viscosity

RHS=(1/k)*sqrt(Cf/2)*log(y.*U1/mu)  + ...
    (1/k)*sqrt(Cf/2)*log(sqrt(Cf/2) + ...
    A * sqrt(Cf/2)



%once Cf is found out, we can calculate tau_w and then U_tau
tau_2=Cf*(1/2)*ro*U1.^2
U_tau= sqrt(tau_w/ro)
