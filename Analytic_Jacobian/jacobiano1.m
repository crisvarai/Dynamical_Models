% JACOBIANO GEMETRICO DEL ROBOT PROMETHEUS
clc
close all
clear all

%----- CINEMATICA DIRECTA DE LOS BRAZOS -------------------

%----- Variables 
syms k q1 q2 q3 q4 q5 q6 d1 a2 a3 d3 a4 d5 d7 w1 w2 w3 w4 w5 w6

%----- Parametros del brazo del Robot Prometheus
%  d1=183.73;
%  a2=58.28;
%  d3=244.4;
%  a3=13.76+5;
%  a4=18.53;
%  d5=78.53+44.5;
%  d7=150;
%  k=((65*pi)/180);

%----- BRAZO IZQUIERDO 

%----- Angulo de giro de cada motor 
%                       %Rango de trabajo
%  q1=0*pi/180;         % -20<q1<180
%  q2=0*pi/180;         % -45<q2<125
%  q3=0*pi/180;         % -45<q3<180
%  q4=0*pi/180;         %  0<q4<135
%  q5=0*pi/180;         %  -180<q5<180
%  q6=0*pi/180;         %  -90<q6<90

%----- Tabla de Denavit-Hartenberg Left 
%           qi   d1   ai     ki
% TL0B=DH(pi/2,  0,   0,    -k);
% TL0H=DH(-pi/2,  0,   0,     0);
% TL01=DH(  q1,  d1,   0, -pi/2);
% TL12=DH(  q2,   0,  a2,  pi/2);
% TL23=DH(  q3,  d3,  a3, -pi/2);
% TL34=DH(  q4,   0,  a4,  pi/2);
% TL45=DH(  q5,  d5,   0, -pi/2);
% TL56=DH(  q6,   0,   0,  pi/2);
% TL67=DH(   0,  d7,   0,     0);

%----- Matriz de orientación-posición del end-efect
%TLee=TL0B*TL0H*TL01*TL12*TL23*TL34*TL45*TL56*TL67

%----- Matrices Homogeneas
TLBH=[cos(k) 0 -sin(k) 0; 0 1 0 0; sin(k) 0 cos(k) 0; 0 0 0 1];
TL01=[cos(q1) 0 -sin(q1) 0; sin(q1) 0 cos(q1) 0; 0 -1 0 d1; 0 0 0 1];
TL12=[cos(q2) 0 sin(q2) a2*cos(q2); sin(q2) 0 -cos(q2) a2*sin(q2); 0 1 0 0; 0 0 0 1];
TL23=[cos(q3) 0 -sin(q3) a3*cos(q3); sin(q3) 0 cos(q3) a3*sin(q3); 0 -1 0 d3; 0 0 0 1];
TL34=[cos(q4) 0 sin(q4) a4*cos(q4); sin(q4) 0 -cos(q4) a4*sin(q4); 0 1 0 0; 0 0 0 1];

%orientación 
TL45=[cos(q5) 0 -sin(q5) 0; sin(q5) 0 cos(q5) 0; 0 -1 0 d5; 0 0 0 1];
TL56=[cos(q6) 0 sin(q6) 0; sin(q1) 0 -cos(q6) 0; 0 1 0 0; 0 0 0 1];
TL67=[1 0 0 0; 0 1 0 0; 0 0 1 d7; 0 0 0 1];

%----- Matriz de orientación-posición del end-efect
TLee=simplify(TLBH*TL01*TL12*TL23*TL34*TL45*TL56*TL67);

%---- JACOBIANO GEOMÉTRICO
% Calculo del Jacobiano de la velocidad lineal Jv
x=TLee(1,4);
y=TLee(2,4);
z=TLee(3,4);

Dx=simplify([diff(x,q1) diff(x,q2) diff(x,q3) diff(x,q4) diff(x,q5) diff(x,q6)]);
Dy=simplify([diff(y,q1) diff(y,q2) diff(y,q3) diff(y,q4) diff(y,q5) diff(y,q6)]);
Dz=simplify([diff(z,q1) diff(z,q2) diff(z,q3) diff(z,q4) diff(z,q5) diff(z,q6)]);

Jv=[Dx; Dy; Dz];

% Calculo del Jacobiano de la velocidad angular Jw
% Matriz de Rotación del end-efect

R=[TLee(1,1) TLee(1,2) TLee(1,3);
   TLee(2,1) TLee(2,2) TLee(2,3);
   TLee(3,1) TLee(3,2) TLee(3,3)];

%--- W=(dR/dt)R^T
%--- Derivadas parciales de la matriz de rotación
%  DR_q1=diff(R,q1);
%  DR_q2=diff(R,q2);
%  DR_q3=diff(R,q3);
%  DR_q4=diff(R,q4);
%  DR_q5=diff(R,q5);
%  DR_q6=diff(R,q6);

%---- Derivada total respecto al tiempo con wi=dqi/dt 
%  DR=DR_q1*w1 + DR_q2*w2+ DR_q3*w3 + DR_q4*w4 + DR_q5*w5 + DR_q6*w6;

%---- Submatriz antisimetrica  
%  W=DR*transpose(R);

%---- Velocidades angulares en las direcciones x, y, z en función 
%---- de las velocides articulares wi i=1,2,3,4,5,6

%  Wx=W(3,2)
%  wy=W(1,3)
%  Wz=W(2,1)

w1=1;
DR_q1=diff(R,q1)*w1;
W1=DR_q1*transpose(R);

w2=1;
DR_q2=diff(R,q2)*w2;
W2=DR_q2*transpose(R);

w3=1;
DR_q3=diff(R,q3)*w3;
W3=DR_q3*transpose(R);

w4=1;
DR_q4=diff(R,q4)*w4;
W4=DR_q4*transpose(R);

w5=1;
DR_q5=diff(R,q5)*w5;
W5=DR_q5*transpose(R);

w6=1;
DR_q6=diff(R,q6)*w6;
W6=DR_q6*transpose(R);

% Wq1=[W1(3,2) W1(1,3) W1(2,1)];
% Wq2=[W2(3,2) W2(1,3) W2(2,1)];
% Wq3=[W3(3,2) W3(1,3) W3(2,1)];
% Wq4=[W4(3,2) W4(1,3) W4(2,1)];
% Wq5=[W5(3,2) W5(1,3) W5(2,1)];
% Wq6=[W6(3,2) W6(1,3) W6(2,1)];

Wx=simplify([W1(3,2) W2(3,2) W3(3,2) W4(3,2) W5(3,2) W6(3,2)]);
Wy=simplify([W1(1,3) W2(1,3) W3(1,3) W4(1,3) W5(1,3) W6(1,3)]);
Wz=simplify([W1(2,1) W2(2,1) W3(2,1) W4(2,1) W5(2,1) W6(2,1)]);
    
Jw=[Wx; Wy; Wz];

%------ Jacobiano Geométrico 

J=[Jv;
   Jw]

 %--- substitución de valores 
 %--- Angulos del espacio articular 
 
   q_1=0;
   q_2=0;
   q_3=0;
   q_4=0;
   q_5=0;
   q_6=0;
    
 % Jacobiano numérico para un punto en el espacio articular (q1,q2,q3,q4,q5,q6)
 J_numerico=subs(J,[d1, a2, d3, a3, a4, d5, d7, k, q1, q2, q3, q4, q5, q6], [183.73, 58.28, 244.4, 18.76, 18.53, 123.03, 150, 1.134, q_1, q_2, q_3, q_4, q_5, q_6]); 
 J_numerico=double(J_numerico)