%syms C_T N c R dC_Tdmu_z sigma C_Lalpha lambda_i0 X_u m rho A omega_0 T

%% Ferrarese2015 Hexa-Copter Values
omega_0 = 461.9230; % , rad/s
rho = 1.2235; % ρ Air Density, kg/m^3
m = 4; % Aircraft Mass, kg
Iyy = 0.044; % Inertia Tensor, kg/m^2
N_rot = 6; % Number of Rotors
N = 2; % Number of Blades of a Rotor
theta_tw = deg2rad(2); % θtw Blade Twist Angle, rad
Cd = 0.003;
I_rotor = 10^-4;
h = -0.3; % height of rotors from the center of gravity of the aircraft: h < 0 if the rotors are placed above the C.G. itself
Ke = 0.005;
gamma = deg2rad(5); % Γ Dihedral Angles of a Rotor, rad

g = 9.81; % Gravity acceleration, m/s^2
Ixx = 0.044; % Inertia Tensor, kg/m^2
Izz = 0.098; % Inertia Tensor, kg/m^2
Ixz = 0; % Inertia Tensor, kg/m^2
R = 0.15; % Rotor Radius, m
theta_c = deg2rad(15); % θc blades pitch of a rotor in the hovering flight condition, rad
C_Lalpha = 5.5; % Clα blade profile lift curve slope, rad^-1
c = 0.04; % Blade Section Chord, m
b = 0.68; % distance between the C.G. of the aircraft and the rotor disk center, m
Ra = 0.01; % , Ω
xi_1 = deg2rad(5); % ξ Tilting Angles of a Rotor, rad
xi_2 = deg2rad(-5); % ξ Tilting Angles of a Rotor, rad
xi_3 = deg2rad(5); % ξ Tilting Angles of a Rotor, rad
xi_4 = deg2rad(-5); % ξ Tilting Angles of a Rotor, rad
xi_5 = deg2rad(5); % ξ Tilting Angles of a Rotor, rad
xi_6 = deg2rad(-5); % ξ Tilting Angles of a Rotor, rad

tau = 1; % Gear Ratio or Time Constant, 

T_0 = (m * g)/N_rot;
v_i0 = 6.1725; % sqrt(T_0/(2 * rho * A)); % induced velocity of each rotor, m/s

delta_1 = deg2rad(0); % δ Azimuthal Angles of a Rotor, rad
delta_2 = deg2rad(60); % δ Azimuthal Angles of a Rotor, rad
delta_3 = deg2rad(120); % δ Azimuthal Angles of a Rotor, rad
delta_4 = deg2rad(180); % δ Azimuthal Angles of a Rotor, rad
delta_5 = deg2rad(240); % δ Azimuthal Angles of a Rotor, rad
delta_6 = deg2rad(300); % δ Azimuthal Angles of a Rotor, rad

eul1 = [gamma xi_1 delta_1];
eul2 = [gamma xi_2 delta_2];
eul3 = [gamma xi_3 delta_3];
eul4 = [gamma xi_4 delta_4];
eul5 = [gamma xi_5 delta_5];
eul6 = [gamma xi_6 delta_6];
T1 = eul2rotm(eul1,'XYZ');
T2 = eul2rotm(eul2,'XYZ');
T3 = eul2rotm(eul3,'XYZ');
T4 = eul2rotm(eul4,'XYZ');
T5 = eul2rotm(eul5,'XYZ');
T6 = eul2rotm(eul6,'XYZ');
%% Equations
sigma = (N * c)/(pi * R); % Rotor Solidity

A = pi * R^2; % Rotor Disk Area, m^2

% v_i0 = sqrt(T_0/(2 * rho * A)) % induced velocity of each rotor, m/s

lambda_i0 = v_i0/(omega_0 * R); % inflow ratio at hover 

dC_Tdmu_z = (2 * sigma * C_Lalpha * lambda_i0)/(16 * lambda_i0 + C_Lalpha * sigma);

dC_pidmu_z = ((-4 * sigma * C_Lalpha)/(16 * lambda_i0 + sigma * C_Lalpha)) * (theta_c/3 - theta_tw/4 - lambda_i0);

eul1 = [gamma xi_1 delta_1];
eul2 = [gamma xi_2 delta_2];
eul3 = [gamma xi_3 delta_3];
eul4 = [gamma xi_4 delta_4];
eul5 = [gamma xi_5 delta_5];
eul6 = [gamma xi_6 delta_6];
% T1 = eul2rotm(eul1,'XYZ');
% T2 = eul2rotm(eul2,'XYZ');
% T3 = eul2rotm(eul3,'XYZ');
% T4 = eul2rotm(eul4,'XYZ');
% T5 = eul2rotm(eul5,'XYZ');
% T6 = eul2rotm(eul6,'XYZ');
T1 = [cos(gamma)*cos(delta_1)                                   cos(gamma)*sin(delta_1)                                     -sin(gamma);
    sin(xi_1)*sin(gamma)*cos(delta_1)-cos(xi_1)*sin(delta_1)    sin(xi_1)*sin(gamma)*sin(delta_1)+cos(xi_1)*cos(delta_1)    sin(xi_1)*cos(gamma);
    cos(xi_1)*sin(gamma)*cos(delta_1)+sin(xi_1)*sin(delta_1)    cos(xi_1)*sin(gamma)*sin(delta_1)-sin(xi_1)*cos(delta_1)    cos(xi_1)*cos(gamma)]; 
% Rotation Matrix for Rotor 1 Orientation

T2 = [cos(gamma)*cos(delta_2)                                   cos(gamma)*sin(delta_2)                                     -sin(gamma);
    sin(xi_2)*sin(gamma)*cos(delta_2)-cos(xi_2)*sin(delta_2)    sin(xi_2)*sin(gamma)*sin(delta_2)+cos(xi_2)*cos(delta_2)    sin(xi_2)*cos(gamma);
    cos(xi_2)*sin(gamma)*cos(delta_2)+sin(xi_2)*sin(delta_2)    cos(xi_2)*sin(gamma)*sin(delta_2)-sin(xi_2)*cos(delta_2)    cos(xi_2)*cos(gamma)]; 
% Rotation Matrix for Rotor 2 Orientation

T3 = [cos(gamma)*cos(delta_3)                                   cos(gamma)*sin(delta_3)                                     -sin(gamma);
    sin(xi_3)*sin(gamma)*cos(delta_3)-cos(xi_3)*sin(delta_3)    sin(xi_3)*sin(gamma)*sin(delta_3)+cos(xi_3)*cos(delta_3)    sin(xi_3)*cos(gamma);
    cos(xi_3)*sin(gamma)*cos(delta_3)+sin(xi_3)*sin(delta_3)    cos(xi_3)*sin(gamma)*sin(delta_3)-sin(xi_3)*cos(delta_3)    cos(xi_3)*cos(gamma)]; 
% Rotation Matrix for Rotor 3 Orientation

T4 = [cos(gamma)*cos(delta_4)                                   cos(gamma)*sin(delta_4)                                     -sin(gamma);
    sin(xi_4)*sin(gamma)*cos(delta_4)-cos(xi_4)*sin(delta_4)    sin(xi_4)*sin(gamma)*sin(delta_4)+cos(xi_4)*cos(delta_4)    sin(xi_4)*cos(gamma);
    cos(xi_4)*sin(gamma)*cos(delta_4)+sin(xi_4)*sin(delta_4)    cos(xi_4)*sin(gamma)*sin(delta_4)-sin(xi_4)*cos(delta_4)    cos(xi_4)*cos(gamma)]; 
% Rotation Matrix for Rotor 4 Orientation

T5 = [cos(gamma)*cos(delta_5)                                   cos(gamma)*sin(delta_5)                                     -sin(gamma);
    sin(xi_5)*sin(gamma)*cos(delta_5)-cos(xi_5)*sin(delta_5)    sin(xi_5)*sin(gamma)*sin(delta_5)+cos(xi_5)*cos(delta_5)    sin(xi_5)*cos(gamma);
    cos(xi_5)*sin(gamma)*cos(delta_5)+sin(xi_5)*sin(delta_5)    cos(xi_5)*sin(gamma)*sin(delta_5)-sin(xi_5)*cos(delta_5)    cos(xi_5)*cos(gamma)]; 
% Rotation Matrix for Rotor 5 Orientation

T6 = [cos(gamma)*cos(delta_6)                                   cos(gamma)*sin(delta_6)                                     -sin(gamma);
    sin(xi_6)*sin(gamma)*cos(delta_6)-cos(xi_6)*sin(delta_6)    sin(xi_6)*sin(gamma)*sin(delta_6)+cos(xi_6)*cos(delta_6)    sin(xi_6)*cos(gamma);
    cos(xi_6)*sin(gamma)*cos(delta_6)+sin(xi_6)*sin(delta_6)    cos(xi_6)*sin(gamma)*sin(delta_6)-sin(xi_6)*cos(delta_6)    cos(xi_6)*cos(gamma)]; 
% Rotation Matrix for Rotor 6 Orientation

eula = [gamma xi_1 0];
eulb = [gamma xi_2 0];
eulc = [gamma xi_3 0];
euld = [gamma xi_4 0];
eule = [gamma xi_5 0];
eulf = [gamma xi_6 0];
T1_tilde = eul2rotm(eula,'XYZ');
T2_tilde = eul2rotm(eulb,'XYZ');
T3_tilde = eul2rotm(eulc,'XYZ');
T4_tilde = eul2rotm(euld,'XYZ');
T5_tilde = eul2rotm(eule,'XYZ');
T6_tilde = eul2rotm(eulf,'XYZ');

T1_tilde = [cos(gamma)              0          -sin(gamma);
           sin(xi_1)*sin(gamma)    cos(xi_1)    sin(xi_1)*cos(gamma);
           cos(xi_1)*sin(gamma)    sin(xi_1)    cos(xi_1)*cos(gamma)]; 
% Rotation Matrix for Rotor 1 Orientation independent of the azimuth angle δ
% i.e set δ = 0

T2_tilde = [cos(gamma)              0          -sin(gamma);
           sin(xi_2)*sin(gamma)    cos(xi_2)    sin(xi_2)*cos(gamma);
           cos(xi_2)*sin(gamma)    sin(xi_2)    cos(xi_2)*cos(gamma)]; 
% Rotation Matrix for Rotor 1 Orientation independent of the azimuth angle δ
% i.e set δ = 0

T3_tilde = [cos(gamma)              0          -sin(gamma);
           sin(xi_3)*sin(gamma)    cos(xi_3)    sin(xi_3)*cos(gamma);
           cos(xi_3)*sin(gamma)    sin(xi_3)    cos(xi_3)*cos(gamma)]; 
% Rotation Matrix for Rotor 1 Orientation independent of the azimuth angle δ
% i.e set δ = 0

T4_tilde = [cos(gamma)              0          -sin(gamma);
           sin(xi_4)*sin(gamma)    cos(xi_4)    sin(xi_4)*cos(gamma);
           cos(xi_4)*sin(gamma)    sin(xi_4)    cos(xi_4)*cos(gamma)]; 
% Rotation Matrix for Rotor 1 Orientation independent of the azimuth angle δ
% i.e set δ = 0

T5_tilde = [cos(gamma)              0          -sin(gamma);
           sin(xi_5)*sin(gamma)    cos(xi_5)    sin(xi_5)*cos(gamma);
           cos(xi_5)*sin(gamma)    sin(xi_5)    cos(xi_5)*cos(gamma)]; 
% Rotation Matrix for Rotor 1 Orientation independent of the azimuth angle δ
% i.e set δ = 0

T6_tilde = [cos(gamma)              0          -sin(gamma);
           sin(xi_6)*sin(gamma)    cos(xi_6)    sin(xi_6)*cos(gamma);
           cos(xi_6)*sin(gamma)    sin(xi_6)    cos(xi_6)*cos(gamma)]; 
% Rotation Matrix for Rotor 1 Orientation independent of the azimuth angle δ
% i.e set δ = 0

% Analytic expressions for stability derivatives

X_u = -1/m * (dC_Tdmu_z * rho * A * omega_0 * R * T1(3,1)^2 + ... 
        dC_Tdmu_z * rho * A * omega_0 * R * T2(3,1)^2 + ...
        dC_Tdmu_z * rho * A * omega_0 * R * T3(3,1)^2 + ...
        dC_Tdmu_z * rho * A * omega_0 * R * T4(3,1)^2 + ...
        dC_Tdmu_z * rho * A * omega_0 * R * T5(3,1)^2 + ...
        dC_Tdmu_z * rho * A * omega_0 * R * T6(3,1)^2);

X_q = -1/m * (dC_Tdmu_z * rho * A * omega_0 * R * (-b * cos(delta_1)) * T1(3,1) * T1(3,3) + ...
    dC_Tdmu_z * rho * A * omega_0 * R * (-b * cos(delta_2)) * T2(3,1) * T2(3,3) + ...
    dC_Tdmu_z * rho * A * omega_0 * R * (-b * cos(delta_3)) * T3(3,1) * T3(3,3) + ...
    dC_Tdmu_z * rho * A * omega_0 * R * (-b * cos(delta_4)) * T4(3,1) * T4(3,3) + ...
    dC_Tdmu_z * rho * A * omega_0 * R * (-b * cos(delta_5)) * T5(3,1) * T5(3,3) + ...
    dC_Tdmu_z * rho * A * omega_0 * R * (-b * cos(delta_6)) * T6(3,1) * T6(3,3)) - ...
     (-1/m) * (dC_Tdmu_z * rho * A * omega_0 * R * -h * T1(3,1)^2 + ...
     dC_Tdmu_z * rho * A * omega_0 * R * -h * T2(3,1)^2 + ...
     dC_Tdmu_z * rho * A * omega_0 * R * -h * T3(3,1)^2 + ...
     dC_Tdmu_z * rho * A * omega_0 * R * -h * T4(3,1)^2 + ...
     dC_Tdmu_z * rho * A * omega_0 * R * -h * T5(3,1)^2 + ...
     dC_Tdmu_z * rho * A * omega_0 * R * -h * T6(3,1)^2);

Y_v = -1/m * (dC_Tdmu_z * rho * A * omega_0 * R * T1(3,2)^2 + ...
    dC_Tdmu_z * rho * A * omega_0 * R * T2(3,2)^2 + ...
    dC_Tdmu_z * rho * A * omega_0 * R * T3(3,2)^2 + ...
    dC_Tdmu_z * rho * A * omega_0 * R * T4(3,2)^2 + ...
    dC_Tdmu_z * rho * A * omega_0 * R * T5(3,2)^2 + ...
    dC_Tdmu_z * rho * A * omega_0 * R * T6(3,2)^2);

Y_p = -1/m * (dC_Tdmu_z * rho * A * omega_0 * R * ((b * sin(delta_1) * T1(3,2) * T1(3,3) + ((-h) * T1(3,2)^2))) + ...
    dC_Tdmu_z * rho * A * omega_0 * R * ((b * sin(delta_2) * T2(3,2) * T2(3,3) + ((-h) * T2(3,2)^2))) + ...
    dC_Tdmu_z * rho * A * omega_0 * R * ((b * sin(delta_3) * T3(3,2) * T3(3,3) + ((-h) * T3(3,2)^2))) + ...
    dC_Tdmu_z * rho * A * omega_0 * R * ((b * sin(delta_4) * T4(3,2) * T4(3,3) + ((-h) * T4(3,2)^2))) + ...
    dC_Tdmu_z * rho * A * omega_0 * R * ((b * sin(delta_5) * T5(3,2) * T5(3,3) + ((-h) * T5(3,2)^2))) + ...
    dC_Tdmu_z * rho * A * omega_0 * R * ((b * sin(delta_6) * T6(3,2) * T6(3,3) + ((-h) * T6(3,2)^2))));

Z_w = -1/m * (dC_Tdmu_z * rho * A * omega_0 * R * T1(3,3)^2 + ...
    dC_Tdmu_z * rho * A * omega_0 * R * T2(3,3)^2 + ...
    dC_Tdmu_z * rho * A * omega_0 * R * T3(3,3)^2 + ...
    dC_Tdmu_z * rho * A * omega_0 * R * T4(3,3)^2 + ...
    dC_Tdmu_z * rho * A * omega_0 * R * T5(3,3)^2 + ...
    dC_Tdmu_z * rho * A * omega_0 * R * T6(3,3)^2);

L_v = -1/Ixx * (dC_Tdmu_z * rho * A * omega_0 * R * T1(3,2) * ((b * sin(delta_1) * T1(3,3) + ((-h) * T1(3,2)))) + ...
    dC_Tdmu_z * rho * A * omega_0 * R * T2(3,2) * ((b * sin(delta_2) * T2(3,3) + ((-h) * T2(3,2)))) + ...
    dC_Tdmu_z * rho * A * omega_0 * R * T3(3,2) * ((b * sin(delta_3) * T3(3,3) + ((-h) * T3(3,2)))) + ...
    dC_Tdmu_z * rho * A * omega_0 * R * T4(3,2) * ((b * sin(delta_4) * T4(3,3) + ((-h) * T4(3,2)))) + ...
    dC_Tdmu_z * rho * A * omega_0 * R * T5(3,2) * ((b * sin(delta_5) * T5(3,3) + ((-h) * T5(3,2)))) + ...
    dC_Tdmu_z * rho * A * omega_0 * R * T6(3,2) * ((b * sin(delta_6) * T6(3,3) + ((-h) * T6(3,2)))));

L_p = -1/Ixx * (dC_Tdmu_z * rho * A * omega_0 * R * ((b * sin(delta_1) * T1(3,3) - h * T1(3,2)))^2 + ...
    dC_Tdmu_z * rho * A * omega_0 * R * ((b * sin(delta_2) * T2(3,3) - h * T2(3,2)))^2 + ...
    dC_Tdmu_z * rho * A * omega_0 * R * ((b * sin(delta_3) * T3(3,3) - h * T3(3,2)))^2 + ...
    dC_Tdmu_z * rho * A * omega_0 * R * ((b * sin(delta_4) * T4(3,3) - h * T4(3,2)))^2 + ...
    dC_Tdmu_z * rho * A * omega_0 * R * ((b * sin(delta_5) * T5(3,3) - h * T5(3,2)))^2 + ...
    dC_Tdmu_z * rho * A * omega_0 * R * ((b * sin(delta_6) * T6(3,3) - h * T6(3,2)))^2);

M_u = -1/Iyy * (dC_Tdmu_z * rho * A * omega_0 * R * T1(3,1) * (h * T1(3,1) + ((-b * cos(delta_1)) * T1(3,3))) + ...
    dC_Tdmu_z * rho * A * omega_0 * R * T2(3,1) * (h * T2(3,1) + ((-b * cos(delta_2)) * T2(3,3))) + ...
    dC_Tdmu_z * rho * A * omega_0 * R * T3(3,1) * (h * T3(3,1) + ((-b * cos(delta_3)) * T3(3,3))) + ...
    dC_Tdmu_z * rho * A * omega_0 * R * T4(3,1) * (h * T4(3,1) + ((-b * cos(delta_4)) * T4(3,3))) + ...
    dC_Tdmu_z * rho * A * omega_0 * R * T5(3,1) * (h * T5(3,1) + ((-b * cos(delta_5)) * T5(3,3))) + ...
    dC_Tdmu_z * rho * A * omega_0 * R * T6(3,1) * (h * T6(3,1) + ((-b * cos(delta_6)) * T6(3,3))));

M_q = -1/Iyy * (dC_Tdmu_z * rho * A * omega_0 * R * ((b * cos(delta_1) * T1(3,3) + ((-h) * T1(3,1))))^2 + ...
    dC_Tdmu_z * rho * A * omega_0 * R * ((b * cos(delta_2) * T2(3,3) + ((-h) * T2(3,1))))^2 + ...
    dC_Tdmu_z * rho * A * omega_0 * R * ((b * cos(delta_3) * T3(3,3) + ((-h) * T3(3,1))))^2 + ...
    dC_Tdmu_z * rho * A * omega_0 * R * ((b * cos(delta_4) * T4(3,3) + ((-h) * T4(3,1))))^2 + ...
    dC_Tdmu_z * rho * A * omega_0 * R * ((b * cos(delta_5) * T5(3,3) + ((-h) * T5(3,1))))^2 + ...
    dC_Tdmu_z * rho * A * omega_0 * R * ((b * cos(delta_6) * T6(3,3) + ((-h) * T6(3,1))))^2);

N_r = -1/Izz * (dC_Tdmu_z * rho * A * omega_0 * R * b^2 * T1_tilde(3,2)^2 + ...
    dC_Tdmu_z * rho * A * omega_0 * R * b^2 * T2_tilde(3,2)^2 + ...
    dC_Tdmu_z * rho * A * omega_0 * R * b^2 * T3_tilde(3,2)^2 + ...
    dC_Tdmu_z * rho * A * omega_0 * R * b^2 * T4_tilde(3,2)^2 + ...
    dC_Tdmu_z * rho * A * omega_0 * R * b^2 * T5_tilde(3,2)^2 + ...
    dC_Tdmu_z * rho * A * omega_0 * R * b^2 * T6_tilde(3,2)^2) + ...
    1/Izz * (dC_pidmu_z * rho * A * omega_0 * +R^2 * b * T1_tilde(3,2) * T1_tilde(3,3) + ...
    dC_pidmu_z * rho * A * omega_0 * -R^2 * b * T2_tilde(3,2) * T2_tilde(3,3) + ...
    dC_pidmu_z * rho * A * omega_0 * +R^2 * b * T3_tilde(3,2) * T3_tilde(3,3) + ...
    dC_pidmu_z * rho * A * omega_0 * -R^2 * b * T4_tilde(3,2) * T4_tilde(3,3) + ...
    dC_pidmu_z * rho * A * omega_0 * +R^2 * b * T5_tilde(3,2) * T5_tilde(3,3) + ...
    dC_pidmu_z * rho * A * omega_0 * -R^2 * b * T6_tilde(3,2) * T6_tilde(3,3));

%% Taken from DeAngelis2018, equation (31)
N_r = -(rho * A * omega_0 * R * b)/Izz * (b * dC_Tdmu_z * ((T1_tilde(3,2)^2)+...
    (T2_tilde(3,2)^2)+...
    (T3_tilde(3,2)^2)+...
    (T4_tilde(3,2)^2)+...
    (T5_tilde(3,2)^2)+...
    (T6_tilde(3,2)^2)) - R * dC_pidmu_z * (((T1_tilde(3,2) * T1_tilde(3,3)) * +1)+...
    ((T2_tilde(3,2) * T2_tilde(3,3)) * -1)+...
    ((T3_tilde(3,2) * T3_tilde(3,3)) * +1)+...
    ((T4_tilde(3,2) * T4_tilde(3,3)) * -1)+...
    ((T5_tilde(3,2) * T5_tilde(3,3)) * +1)+...
    ((T6_tilde(3,2) * T6_tilde(3,3)) * -1)));

% GIVES SAME RESULT

X_theta = -g;

Y_phi = g;

% X_u = calculated above
% X_q = calculated above
% Y_v = calculated above
% Y_p = calculated above
% Z_w = calculated above
% L_v = calculated above
% L_p = calculated above
% M_u = calculated above
% M_q = calculated above
% N_r = calculated above

X_w = 0; % Xw is null because the variation of C_T is null. 
% The rotors drag effects due to yaw rate variations are nullified by the
% symmetric displacement of the rotors. (Ferrarese2015)
Z_u = 0; % Along the zB axis the perturbations of forces have zero contributions from u because
% deltaCT/deltamu = 0, equation 5.11, (Ferrarese2015)
Z_q = 0; % Along the zB axis, For symmetry in rotors displacement the effects of q is null. (Ferrarese2015) 
M_w = 0; % For reasons of symmetry, other effects are negligible. (Ferrarese2015)
Y_r = 0; % The rotors drag e?ects due to yaw rate variations are nullifed by the symmetric displacement of the rotors. (Ferrarese2015)
L_r = 0; % For reasons of symmetry, other effects are negligible. (Ferrarese2015)
N_v = 0; % The other effects can be supposed negligible for symmetry of rotors displacement 
% and for the alternation of their verses of rotation. (Ferrarese2015)
N_p = 0; % The other effects can be supposed negligible for symmetry of rotors displacement 
% and for the alternation of their verses of rotation. (Ferrarese2015)


AF = [0     0       0   0       0   0   1   0   0;
    0       0       0   0       0   0   0   1   0;
    0       0       0   0       0   0   0   0   1;
    0       X_theta 0   X_u     0   0   0   X_q 0;
    Y_phi   0       0   0       Y_v 0   Y_p 0   0;
    0       0       0   0       0   Z_w 0   0   0;
    0       0       0   0       L_v 0   L_p 0   0;
    0       0       0   M_u     0   0   0   M_q 0;
    0       0       0   0       0   0   0   0   N_r]
    
%% Longitudinal-Heave
% xlongF = [delta_u; delta_w; delta_q; delta_theta]

AlongF = [X_u/m     X_w/m       X_q/m       -g;
      Z_u/m         Z_w/m       Z_q/m       0;
      M_u/Iyy        M_w/Iyy      M_q/Iyy      0;
      0             0           1           0]

[V,D] = eig(AlongF)  
  
AlongF1 = [X_u      X_w         X_q         -g;
      Z_u           Z_w         Z_q         0;
      M_u           M_w         M_q         0;
      0             0           1           0]
  
[V,D] = eig(AlongF1)
 %% Lateral-Directional
% xlatF = [delta_v; delta_p; delta_r; delta_phi]

AlatF = [Y_v/m                          Y_p/m                               Y_r/m                               g;
      (Izz*L_v+Ixz*N_v)/(Ixx*Izz-Ixz^2)    (Izz*L_p+Ixz*N_p)/(Ixx*Izz-Ixz^2)      (Izz*L_r+Ixz*N_r)/(Ixx*Izz-Ixz^2)      0;
      (Ixz*L_v+Ixx*N_v)/(Ixx*Izz-Ixz^2)    (Ixz*L_p+Ixx*N_p)/(Ixx*Izz-Ixz^2)      (Ixz*L_r+Ixx*N_r)/(Ixx*Izz-Ixz^2)      0;
      0                                 1                                   0                                   0]

[V,D] = eig(AlatF)    

AlatF1 = [Y_v    Y_p     Y_r     g;
      L_v    L_p      L_r      0;
      N_v    N_p      N_r      0;
      0     1       0           0] 

[V,D] = eig(AlatF1)    