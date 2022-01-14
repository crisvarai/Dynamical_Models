%% Ley de Control

clc;        % Limpiar Command Window
clear all   % Limpiar Workspace
close all   % Cerrar todas las ventanas de simulación
disp('*** SIMULACION DE LA LEY DE CONTROL  ***');

%% Declaración de variables

syms Ts % Tiempo de muestreo
syms xr yr zr nr br1 br2 br3 % Pose del robot
syms xl yl zl nl bl1 bl2 bl3 % Pose del operador
syms ypxr ypyr ypzr ypnr ypbr1 ypbr2 ypbr3 % Dinámica del filtro
syms yxr yyr yzr ynr ybr1 ybr2 ybr3 % Integral del filtro
syms kyi % Ganancia del filtro

syms kr1 kr2 kr3 kr4 kr5 kr6 % Ganancia proporcional
syms dr1 dr2 dr3 dr4 dr5 dr6 % Ganancia derativa

% Cálculo de Velocidades
% ypxr = xr-kyi*yxr;
% yxr = yxr + Ts*ypxr;
% 
% ypyr = yr-kyi*yyr;
% yyr = yyr + Ts*ypyr;
% 
% ypzr = zr-kyi*yzr;
% yzr = yzr + Ts*ypzr;
% 
% ypnr = nr-kyi*ynr;
% ynr = ynr + Ts*ypnr;
% 
% ypbr1 = br1-kyi*ybr1;
% ybr1 = ybr1 + Ts*ypbr1;
% 
% ypbr2 = br2-kyi*ybr2;
% ybr2 = ybr2 + Ts*ypbr2;
% 
% ypbr3 = br3-kyi*ybr3;
% ybr3 = ybr3 + Ts*ypbr3;

%% Cálculo del termino A1

% Ganancia proporcional
Kr = [kr1; kr2; kr3; kr4; kr5; kr6];

% Cálculo del error
pr = [xr; yr; zr];
pl = [xl; yl; zl];
br = [br1; br2; br3];
bl = [bl1; bl2; bl3];
S_br = [0 -br3 br2; br3 0 -br1; -br2 br1 0];

er = [pr - pl; nl*br - nr*bl - S_br*bl];

% Cálculo de A1
A1 = -[Kr(1)*er(1); Kr(2)*er(2); Kr(3)*er(3); Kr(4)*er(4); Kr(5)*er(5); Kr(6)*er(6)];

%% Cálculo del termino A2

% Ganancia proporcional
Dr = [dr1; dr2; dr3; dr4; dr5; dr6];

% Cálculo de U_Ei
U_Ei = [-transpose(br); nr*eye(3)-S_br];
UT_Ei = transpose(U_Ei);

% Cálculo de Psi
psi_r = [eye(3) zeros(3,4); zeros(3) 0.5*(UT_Ei)];

% Cálculo de yp_i (velocidades)
yp_i = [ypxr; ypyr; ypzr; ypnr; ypbr1; ypbr2; ypbr3];

% Cálculo de A1
A3_inc = [yp_i(1)*psi_r(:,1) + yp_i(2)*psi_r(:,2) + yp_i(3)*psi_r(:,3) + yp_i(4)*psi_r(:,4) + yp_i(5)*psi_r(:,5) + yp_i(6)*psi_r(:,6) + yp_i(7)*psi_r(:,7)];
A3 = -[Dr(1)*A3_inc(1); Dr(2)*A3_inc(2); Dr(3)*A3_inc(3); Dr(4)*A3_inc(4); Dr(5)*A3_inc(5); Dr(6)*A3_inc(6)];

%% Cálculo del vector de Fuerzas

f = A1+A3