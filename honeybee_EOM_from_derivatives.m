%% Symbolic Expressions

syms m_n_j m_n Ix_n Iy_n_j Iy_n Iz_n Ixz_n g_n_j g_n...
    X_u_n X_w_n X_q_n Z_u_n Z_w_n Z_q_n M_u_n M_w_n M_q_n...
    Y_v_n L_v_n N_v_n Y_p_n L_p_n N_p_n Y_r_n L_r_n N_r_n...
    delta_u_n delta_w_n delta_q_n delta_theta_n...
    delta_v_n delta_p_n delta_r_n delta_phi_n

%% Jiang2007 - Longitudinal-Heave
xlongJ_n = [delta_u_n; delta_w_n; delta_q_n; delta_theta_n]

AlongJ_n = [X_u_n/m_n_j     X_w_n/m_n_j       X_q_n/m_n_j       -g_n_j;
      Z_u_n/m_n_j           Z_w_n/m_n_j       Z_q_n/m_n_j       0;
      M_u_n/Iy_n_j          M_w_n/Iy_n_j     M_q_n/Iy_n_j       0;
      0                     0               1                   0]
  
AlongJ_n*xlongJ_n  
 %% Xu2014 - Lateral-Directional
xlatX_n = [delta_v_n; delta_p_n; delta_r_n; delta_phi_n]

AlatX_n = [Y_v_n/m_n                                      Y_p_n/m_n                                       Y_r_n/m_n                                       g_n;
      (Iz_n*L_v_n+Ixz_n*N_v_n)/(Ix_n*Iz_n-Ixz_n^2)      (Iz_n*L_p_n+Ixz_n*N_p_n)/(Ix_n*Iz_n-Ixz_n^2)    (Iz_n*L_r_n+Ixz_n*N_r_n)/(Ix_n*Iz_n-Ixz_n^2)    0;
      (Ixz_n*L_v_n+Ix_n*N_v_n)/(Ix_n*Iz_n-Ixz_n^2)      (Ixz_n*L_p_n+Ix_n*N_p_n)/(Ix_n*Iz_n-Ixz_n^2)    (Ixz_n*L_r_n+Ix_n*N_r_n)/(Ix_n*Iz_n-Ixz_n^2)    0;
      0                                                 1                                               0                                               0] 

AlatX_n*xlatX_n
