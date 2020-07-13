syms u v w p q r phi theta psi...
    X_u X_v X_w X_p X_q W_0 X_r V_0 X_theta...
    g theta_0 Y_u Y_v Y_w Y_p W_0 Y_q Y_r U_0 Y_phi g...
    Z_u Z_v Z_w Z_p V_0 Z_q Z_r... 
    L_u L_v L_w L_p L_q L_r L_phi... 
    M_u M_v M_w M_p M_q M_r M_theta... 
    N_u N_v N_w N_p N_q N_r...
    m I_x I_y I_z I_xz

m = 101.9/(1000*1000); % total mass of insect, kg
mw = 0.50; % percentage of the mass of wing pair,; %
R = 9.8/1000; % wing length, m
c = 2.91/1000; % mean chord length of wing, m
% c = 3.08/1000; % mean chord length of wing, m Jiang2007
Ix = 0.943*10^-9; % moment of inertia of the insect (body and wings) about x- axis, kg·m^2
Iy = 1.621*10^-9; % moment of inertia of the insect (body and wings) about y- axis, kg·m^2
Iz = 0.862*10^-9; % moment of inertia of the insect (body and wings) about x- and z- axes, kg·m^2
Ixz = -0.716*10^-9; % product of inertia of the insect (body and wings) about x- and z- axes, kg·m^2
phi = deg2rad(131);% Φ
rho = 1.23;% ρ, kg/m^3
U = 2 * phi * n * r^2; % mean flapping velocity
tw = 1/n; % , the period of the wingbeat cycle
St = ; % the area of two wings
G = 9.81; % m/s^2
m_n = m/(0.5 * rho * U * St * tw);
Ix_n = Ix/(0.5 * rho * U^2 * St * c * tw^2);
Iy_n = Iy/(0.5 * rho * U^2 * St * c * tw^2); 
Iz_n = Iz/(0.5 * rho * U^2 * St * c * tw^2); 
Ixz_n = Ixz/(0.5 * rho * U^2 * St * c * tw^2); 
g_n = g * tw/U;

delta_u_n = delta_u/U;
delta_w_n = delta_w/U;
delta_q_n = delta_q * tw;
delta_v_n = delta_v/U; 
delta_p_n = delta_p * tw;
delta_r_n = delta_r * tw;


X_u_n
X_w_n
X_q_n
Z_u_n
Z_w_n
Z_q_n
M_u_n
M_w_n
M_q_n
Y_v_n = -0.99;
L_v_n = -0.63;
N_v_n = -0.07;
Y_p_n = -0.10;
L_p_n = -1.09;
N_p_n = -0.06;
Y_r_n = 0;
L_r_n = 0.01;
N_r_n = -1.09;

X_n = X/(0.5 * rho * U^2 * St); 
Z_n = Z/(0.5 * rho * U^2 * St); 
M_n = M/(0.5 * rho * U^2 * St * c);
Y_n = Y/(0.5 * rho * U^2 * St); 
L_n = L/(0.5 * rho * U^2 * St * c); 
N_n = N/(0.5 * rho * U^2 * St * c); 
t_n = t/tw;

X = X_n*(0.5 * rho * U^2 * St); 
Z = Z_n*(0.5 * rho * U^2 * St); 
M = M_n*(0.5 * rho * U^2 * St * c); 
Y = Y_n*(0.5 * rho * U^2 * St); 
L = L_n*(0.5 * rho * U^2 * St * c); 
N = N_n*(0.5 * rho * U^2 * St * c); 
%% Jiang2007
xJ = [u; w; q; theta]


AlongJ_n = [X_u_n/m_n     X_w_n/m_n       X_q_n/m_n       -g_n;
      Z_u_n/m_n         Z_w_n/m_n       Z_q_n/m_n       0;
      M_u_n/I_y_n       M_w_n/I_y_n     M_q_n/I_y_n     0;
      0                 0               1               0]

AlongJ_n*xJ

 %% Xu2014
xX = [v; p; r; phi]


AlatX_n = [Y_v_n/m_n                                      Y_p_n/m_n                                       Y_r_n/m_n                                       g_n;
      (Iz_n*L_v_n+Ixz_n*N_v_n)/(Ix_n*Iz_n-Ixz_n^2)      (Iz_n*L_p_n+Ixz_n*N_p_n)/(Ix_n*Iz_n-Ixz_n^2)    (Iz_n*L_r_n+Ixz_n*N_r_n)/(Ix_n*Iz_n-Ixz_n^2)    0;
      (Ixz_n*L_v_n+Ix_n*N_v_n)/(Ix_n*Iz_n-Ixz_n^2)      (Ixz_n*L_p_n+Ix_n*N_p_n)/(Ix_n*Iz_n-Ixz_n^2)    (Ixz_n*L_r_n+Ix_n*N_r_n)/(Ix_n*Iz_n-Ixz_n^2)    0;
      0                                                 1                                               0                                               0] 
  
AlatX_n*xX

