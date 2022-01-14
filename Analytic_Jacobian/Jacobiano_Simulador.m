
% run startup_rvc   Correr 1 vez para cargar la librería

clc;        % Limpiar Command Window
clear all   % Limpiar Workspace
close all   % Cerrar todas las ventanas de simulación
disp('*** JACOBIANO ANALITICO - BRAZO IZQUIERDO  ***');

%% ---------------------- Información del Script --------------------------

% PROGRAMA: Que efectúa la simulación del Jacobiano Analítico de un brazo 
%           (izquierdo) robótico antropomórfico de 6DOF 
% OBJETIVO: Validar el Jacobiano Analítico calculado teóricamente
% FECHA:    9 de Julio 2021
% DISEÑO:   Ing. Cristian Vallejo

%% ------------ Parámetros DH del Brazo Izquierdo del Robot ---------------

d1 = 183.73/1000;
a2 = 58.28/1000;
d3 = 244.4/1000;
a3 = 18.76/1000;
a4 = 18.53/1000;
d5 = 123.03/1000;
d7 = 150/1000;

alpha = 65;
alpha = alpha*pi/180;

offset = 90;
offset = offset*pi/180;

theta = [0 0 0 0 0 0];
% theta = []
% count = 0;
% while 1
%     fprintf('Ingresa theta %d: ', count+1);
%     theta(count+1) = input('');
%     count = count+1;
%     if count == 6; break; end
% end
theta = theta*pi/180;

%% -------- Implementación de los parámetros DH en el Serial-Link ---------

L(1) = Link('d', 0,  'a', 0,  'alpha', -alpha, 'offset',  offset);
L(2) = Link('d', 0,  'a', 0,  'alpha',  0,     'offset', -offset);
L(3) = Link('d', d1, 'a', 0,  'alpha', -pi/2,  'offset',  theta(1));
L(4) = Link('d', 0,  'a', a2, 'alpha',  pi/2,  'offset',  theta(2));
L(5) = Link('d', d3, 'a', a3, 'alpha', -pi/2,  'offset',  theta(3));
L(6) = Link('d', 0,  'a', a4, 'alpha',  pi/2,  'offset',  theta(4));
L(7) = Link('d', d5, 'a', 0,  'alpha', -pi/2,  'offset',  theta(5));
L(8) = Link('d', 0,  'a', 0,  'alpha',  pi/2,  'offset',  theta(6));
L(9) = Link('d', d7, 'a', 0,  'alpha',  0,     'offset',  0);

LeftArm = SerialLink(L);

q0 = [0 0 0 0 0 0 0 0 0];
q0 = q0*pi/180;

qf = [0 0 0 0 0 0 0 0 0];
qf = qf*pi/180;

t = 0:0.15:3;
Q = jtraj(q0,qf,t);

% Jf = [];
% for i=1:1:21
FK = fkine(LeftArm,Q)
J = jacob0(LeftArm,Q(21,1:9));
J = J(1:6,3:8)
% Jf = [Jf;J];
% end

%% --------------------- Simulación de Brazo Derecho ----------------------

LeftArm.name = 'Prometheus Left Arm';
LeftArm.plot([ 0 0 0 0 0 0 0 0 0 ]);
LeftArm.teach();
