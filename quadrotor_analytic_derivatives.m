%syms C_T N c R dC_Tdmu_z sigma C_Lalpha lambda X_u m rho A omega_0 T

%% Ferrarese2015 Hexa-Copter Values
C_Lalpha = 5.5 
m = 4
rho = 1.2235
omega_0 = 461.9230
R = 0.15
T = 1
N = 6
c = 0.04
C_T = 
%% Equations
lambda = sqrt(C_T/2)

sigma = (N*c)/(pi*R)

A = pi*R^2

dC_Tdmu_z = (2*sigma*C_Lalpha*lambda)/(16*lambda+C_Lalpha*sigma)

X_u = 4*(-1/m*dC_Tdmu_z*rho*A*omega_0*R*T^2)

X_u