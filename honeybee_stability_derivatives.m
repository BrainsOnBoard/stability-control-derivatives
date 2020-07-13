%% Parameters
m = 101.9/(1000*1000); % total mass of insect, kg
mw = 0.50; % percentage of the mass of wing pair,; %
R = 9.8/1000; % wing length, m
c_j = 3.08/1000; % mean chord length of wing, m Jiang2007
c = 2.91/1000; % mean chord length of wing, m
r2 = 0.544*R;
Ix = 0.943*10^-9; % moment of inertia of the insect (body and wings) about x-axis, kg·m^2
Iy_j = 1.621*10^-9; % moment of inertia of the insect (body and wings) about y-axis, kg·m^2 Jiang2007
Iy = 1.59*10^-9; % moment of inertia of the insect (body and wings) about y-axis, kg·m^2
Iz = 0.862*10^-9; % moment of inertia of the insect (body and wings) about zaxis, kg·m^2
Ixz = -0.716*10^-9; % product of inertia of the insect (body and wings) about x- and z- axes, kg·m^2
phi = deg2rad(131);% Φ
rho_j = 1.23;% ρ, kg/m^3 Jiang2007
rho = 1.25;% ρ, kg/m^3
n = 197; % Hz
U_j = 5.00; % mean flapping velocity, m/s Jiang2007
U = 2 * phi * n * r2; % mean flapping velocity, m/s
tw = 1/n; % , the period of the wingbeat cycle, s
St_j = 2*R*c_j; % the area of two wings, m^2 Jiang2007 
St = 2*R*c; % the area of two wings, m^2
g = 9.81; % m/s^2

%% Non Dimensional Variables and Derivatives
m_n_j = 108.23;
m_n = m/(0.5 * rho * U * St * tw);
Ix_n = Ix/(0.5 * rho * U^2 * St * c * tw^2);
Iy_n_j = 21.62;
Iy_n = Iy/(0.5 * rho * U^2 * St * c * tw^2); 
Iz_n = Iz/(0.5 * rho * U^2 * St * c * tw^2); 
Ixz_n = Ixz/(0.5 * rho * U^2 * St * c * tw^2); 
g_n_j = 0.01;
g_n = g * tw/U;

% delta_u_n = delta_u/U;
% delta_w_n = delta_w/U;
%delta_q_n = delta_q * tw;
% delta_v_n = delta_v/U; 
% delta_p_n = delta_p * tw;
% delta_r_n = delta_r * tw;


X_u_n = -0.0057*m_n_j;
X_w_n = 0.0002*m_n_j;
X_q_n = -0.0002*m_n_j;
Z_u_n = -0.0007*m_n_j;
Z_w_n = -0.0100*m_n_j;
Z_q_n = 0.0001*m_n_j;
M_u_n = 0.1295*Iy_n_j;
M_w_n = 0.0035*Iy_n_j;
M_q_n = -0.0675*Iy_n_j;
Y_v_n = -0.99;
L_v_n = -0.63;
N_v_n = -0.07;
Y_p_n = -0.10;
L_p_n = -1.09;
N_p_n = -0.06;
Y_r_n = 0;
L_r_n = 0.01;
N_r_n = -1.09;

% X_n = X/(0.5 * rho * U^2 * St); 
% Z_n = Z/(0.5 * rho * U^2 * St); 
% M_n = M/(0.5 * rho * U^2 * St * c);
% Y_n = Y/(0.5 * rho * U^2 * St); 
% L_n = L/(0.5 * rho * U^2 * St * c); 
% N_n = N/(0.5 * rho * U^2 * St * c); 
% t_n = t/tw;

%% Computed Pure Stability Derivatives
X_u = X_u_n*(0.5 * rho_j * U_j^2 * St_j); 
X_w = X_w_n*(0.5 * rho_j * U_j^2 * St_j); 
X_q = X_q_n*(0.5 * rho_j * U_j^2 * St_j); 
Z_u = Z_u_n*(0.5 * rho_j * U_j^2 * St_j); 
Z_w = Z_w_n*(0.5 * rho_j * U_j^2 * St_j); 
Z_q = Z_q_n*(0.5 * rho_j * U_j^2 * St_j); 
M_u = M_u_n*(0.5 * rho_j * U_j^2 * St_j * c_j); 
M_w = M_w_n*(0.5 * rho_j * U_j^2 * St_j * c_j); 
M_q = M_q_n*(0.5 * rho_j * U_j^2 * St_j * c_j); 
Y_v = Y_v_n*(0.5 * rho * U^2 * St); 
Y_p = Y_p_n*(0.5 * rho * U^2 * St); 
Y_r = Y_r_n*(0.5 * rho * U^2 * St); 
L_v = L_v_n*(0.5 * rho * U^2 * St * c);
L_p = L_p_n*(0.5 * rho * U^2 * St * c);
L_r = L_r_n*(0.5 * rho * U^2 * St * c);
N_v = N_v_n*(0.5 * rho * U^2 * St * c);
N_p = N_p_n*(0.5 * rho * U^2 * St * c);
N_r = N_r_n*(0.5 * rho * U^2 * St * c);

%% Jiang2007 - Longitudinal-Heave
% xlongJ_n = [delta_u_n; delta_w_n; delta_q_n; delta_theta_n]

AlongJ_n = [X_u_n/m_n_j     X_w_n/m_n_j       X_q_n/m_n_j       -g_n_j;
      Z_u_n/m_n_j           Z_w_n/m_n_j       Z_q_n/m_n_j       0;
      M_u_n/Iy_n_j          M_w_n/Iy_n_j     M_q_n/Iy_n_j       0;
      0                     0               1                   0]
  
% xlongJ = [delta_u; delta_w; delta_q; delta_theta]

[V,D] = eig(AlongJ_n)

AlongJ = [X_u/m     X_w/m       X_q/m       -g;
      Z_u/m         Z_w/m       Z_q/m       0;
      M_u/Iy        M_w/Iy      M_q/Iy      0;
      0             0           1           0]
  
[V,D] = eig(AlongJ)
 %% Xu2014 - Lateral-Directional
% xlatX_n = [delta_v_n; delta_p_n; delta_r_n; delta_phi_n]

AlatX_n = [Y_v_n/m_n                                      Y_p_n/m_n                                       Y_r_n/m_n                                       g_n;
      (Iz_n*L_v_n+Ixz_n*N_v_n)/(Ix_n*Iz_n-Ixz_n^2)      (Iz_n*L_p_n+Ixz_n*N_p_n)/(Ix_n*Iz_n-Ixz_n^2)    (Iz_n*L_r_n+Ixz_n*N_r_n)/(Ix_n*Iz_n-Ixz_n^2)    0;
      (Ixz_n*L_v_n+Ix_n*N_v_n)/(Ix_n*Iz_n-Ixz_n^2)      (Ixz_n*L_p_n+Ix_n*N_p_n)/(Ix_n*Iz_n-Ixz_n^2)    (Ixz_n*L_r_n+Ix_n*N_r_n)/(Ix_n*Iz_n-Ixz_n^2)    0;
      0                                                 1                                               0                                               0] 

  
[V,D] = eig(AlatX_n)

% xlatX = [delta_v; delta_p; delta_r; delta_phi]

AlatX = [Y_v/m                          Y_p/m                               Y_r/m                               g;
      (Iz*L_v+Ixz*N_v)/(Ix*Iz-Ixz^2)    (Iz*L_p+Ixz*N_p)/(Ix*Iz-Ixz^2)      (Iz*L_r+Ixz*N_r)/(Ix*Iz-Ixz^2)      0;
      (Ixz*L_v+Ix*N_v)/(Ix*Iz-Ixz^2)    (Ixz*L_p+Ix*N_p)/(Ix*Iz-Ixz^2)      (Ixz*L_r+Ix*N_r)/(Ix*Iz-Ixz^2)      0;
      0                                 1                                   0                                   0] 

[V,D] = eig(AlatX)