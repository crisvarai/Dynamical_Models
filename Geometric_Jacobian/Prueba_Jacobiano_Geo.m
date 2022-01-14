
clc;        % Limpiar Command Window
clear all   % Limpiar Workspace
close all   % Cerrar todas las ventanas de simulación
disp('*** PRUEBA - JACOBIANO GEOMETRICO  ***');

%% -------------- Parámetros DH del Brazo Derecho del Robot ---------------

d1 = 183.73;
a2 = 58.28;
d3 = 244.4;
a3 = 18.76;
a4 = 18.53;
d5 = 123.03;
d7 = 150;

alpha = 65;
alpha = alpha*pi/180;

offset = 90;
offset = offset*pi/180;

%% Cálculo de Transformadas

t1 = 0;
t2 = 0;
t3 = 0;
t4 = 0;
t5 = 0;
t6 = 0;

% Inserting D-H convention parameters
% WAM
A1 = Trans(0,-pi/2,0,t1);
A2 = Trans(0.5518,0,0,t2);
A3 = Trans(-0.45,pi/2,0,t3);
A4 = Trans(0,-pi/2,0.3,t4); 
A5 = Trans(0,pi/2,0,t5);
A6 = Trans(0,0,0,t6);

% Creating Transfer matrices
T2 = A1*A2;
T3 = A1*A2*A3; 
T4 = A1*A2*A3*A4; 
T5 = A1*A2*A3*A4*A5; 
T6 = A1*A2*A3*A4*A5*A6;

% Creating zi 
z0= [0;0;1];  
z1= A1(1:3,3);  
z2= T2(1:3,3);  
z3= T3(1:3,3);  
z4= T4(1:3,3);  
z5= T5(1:3,3);

% Creating pi 
p0=[0;0;0]; 
p1=A1(1:3,4); 
p2=T2(1:3,4); 
p3=T3(1:3,4); 
p4=T4(1:3,4); 
p5=T5(1:3,4); 
P=T6(1:3,4);

% Jacobian matrix Computation
J = ([cross(z0,P-p0),cross(z1,P-p1),cross(z2,P-p2),cross(z3,P-p3),cross(z4,P-p4),cross(z5,P-p5);
    z0,z1,z2,z3,z4,z5]) 

% (3*3) blocks Jacobians 
J11=J(1:3,1:3);
J22=J(4:6,4:6);

%%%%%%% Trans.m %%%%%%%
function [ T ] = Trans(a,b,c,d)
% D-H Homogeneous Transformation Matrix (a alpha d theta)
T = [
cos(d) -sin(d)*round(cos(b)) sin(d)*sin(b) a*cos(d);
sin(d) cos(d)*round(cos(b)) -cos(d)*sin(b) a*sin(d);
0 sin(b) round(cos(b)) c;
0 0 0 1
    ];
end