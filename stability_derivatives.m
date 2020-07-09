
syms u v w p q r phi theta psi...
    X_u X_v X_w X_p X_q W_0 X_r V_0 X_theta...
    g theta_0 Y_u Y_v Y_w Y_p W_0 Y_q Y_r U_0 Y_phi g...
    Z_u Z_v Z_w Z_p V_0 Z_q Z_r... 
    L_u L_v L_w L_p L_q L_r L_phi... 
    M_u M_v M_w M_p M_q M_r M_theta... 
    N_u N_v N_w N_p N_q N_r...
    m I_x I_y I_z I_xz

x = [u; v; w; p; q; r; phi; theta; psi];
%% Niemiec2016
AN = [X_u    0   0   0   0   0       0       X_theta 0;
    0       Y_v 0   0   0   0       Y_phi   0       0;
    0       0   Z_w 0   0   0       0       0       0;
    0       L_v 0   L_p 0   0       L_phi   0       0;
    M_u     0   0   0   M_q 0       0       M_theta 0;
    0       0   0   0   0   N_r     0       0       0;
    0       0   0   1   0   0       0       0       0;
    0       0   0   0   1   0       0       0       0;
    0       0   0   0   0   1       0       0       0];


%% Gong2019 identified the parameters
AG = [X_u    X_v     X_w X_p     X_q-W_0 X_r+V_0         0               -g*cos(theta_0) 0;
    Y_u     Y_v     Y_w Y_p+W_0 Y_q     Y_r-U_0         g*cos(theta_0)  0               0;
    Z_u     Z_v     Z_w Z_p-V_0 Z_q+U_0 Z_r             0               -g*sin(theta_0) 0;
    L_u     L_v     L_w L_p     L_q     L_r             0               0               0;
    M_u     M_v     M_w M_p     M_q     M_r             0               0               0;
    N_u     N_v     N_w N_p     N_q     N_r             0               0               0;
    0       0       0   1       0       tan(theta_0)    0               0               0;
    0       0       0   0       1       0               0               0               0;
    0       0       0   0       0       sec(theta_0)    0               0               0];

%% Ivler2019 similar to Gong2019 but splits longitudinal and lateral
AI = [X_u    0       X_w 0       X_q-W_0 0               0               -g*cos(theta_0) 0;
    0       Y_v     0   0+W_0   0       0-U_0           g               0               0;
    Z_u     0       Z_w 0       Z_q+U_0 0               0               0               0;
    0       L_v     0   L_p     0       L_r             0               0               0;
    M_u     0       M_w 0       M_q     0               0               0               0;
    0       N_v     0   N_p     0       N_r             0               0               0;
    0       0       0   1       0       0               0               0               0;
    0       0       0   0       1       0               0               0               0;
    0       0       0   0       0       1               0               0               0];

% Stability derivatives estimated using RMAC

%% Bristeau determined that without coupling aerodynamic effects to the body, stability derivatives in the matrices above vanish to zero.

%% Craig2019a estimated longitudinal stability derivatives using CIFER.

%% Wei2015, Wei2015a estimated stability derivatives using CIFER.
AW = [X_u    0   0   0   X_q 0       0       -g      0;
    0       Y_v 0   Y_p 0   0       g       0       0;
    0       0   Z_w 0   0   0       0       0       0;
    0       L_v 0   L_p 0   0       0       0       0;
    M_u     0   0   0   M_q 0       0       0       0;
    0       0   0   0   0   N_r     0       0       0;
    0       0   0   1   0   0       0       0       0;
    0       0   0   0   1   0       0       0       0;
    0       0   0   0   0   1       0       0       0];

%% Saetti2018 estimated stability derivatives using CIFER.
AS = [X_u    0   0   0   X_q 0       0       -g      0;
    0       Y_v 0   Y_p 0   0       g       0       0;
    0       0   Z_w 0   0   0       0       0       0;
    0       L_v 0   L_p 0   0       0       0       0;
    M_u     0   0   0   M_q 0       0       0       0;
    0       0   0   0   0   N_r     0       0       0;
    0       0   0   1   0   0       0       0       0;
    0       0   0   0   1   0       0       0       0;
    0       0   0   0   0   1       0       0       0];

% Analytical expressions for stability derivatives for helicopters available from Padfield, Prouty and Wayne Johnson

%% Ferrarese2015
AF = [X_u    0   0   0   0   0       0       X_theta 0;
    0       Y_v 0   0   0   0       Y_phi   0       0;
    0       0   Z_w 0   0   0       0       0       0;
    0       L_v 0   L_p 0   0       0       0       0;
    M_u     0   0   0   M_q 0       0       0       0;
    0       0   0   0   0   N_r     0       0       0;
    0       0   0   1   0   0       0       0       0;
    0       0   0   0   1   0       0       0       0;
    0       0   0   0   0   1       0       0       0];

%% Jiang2007
xj = [u; w; q; theta]


AJ = [X_u/m     X_w/m       X_q/m       -g;
      Z_u/m     Z_w/m       Z_q/m       0;
      M_u/I_y   M_w/I_y     M_q/I_y     0;
      0         0           1           0]
AJ*xj

 %% Xu2014
xx = [v; p; r; phi]


AX = [Y_v/m             Y_p/m               Y_r/m               g;
      N_v               N_p                 N_r                 0;
      N_v               N_p                 N_r                 0;
      0                 1                   1                   0] 
  
AX*xx

