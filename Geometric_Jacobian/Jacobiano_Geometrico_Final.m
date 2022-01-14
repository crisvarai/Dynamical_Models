clear all 
close all 
clc

%%%En este scrip se calcula el Jacobiano Geometrico (Producto cruz) del brazo izquierdo del
%%%robot Prometheus, 

syms q1 q2 q3 q4 q5 q6 %Definimos las variables articulares (simbolicos)
syms k d1 a2 a3 d3 a4 d5 d7 %Definimos las variables de parametros fisicos (simbolicos)

% %%
% %----- Parametros del brazo del Robot Prometheus
%  d1=183.73/1000;
%  a2=58.28/1000;
%  d3=244.4/1000; 
%  a3=13.76+5/1000;
%  a4=18.53/1000;
%  d5=78.53+44.5/1000;
%  d7=150/1000;
%  k=((65*pi)/180);
%  %----- Angulo de giro de cada motor 
% %                       %Rango de trabajo
%   q1=90*pi/180;         % -20<q1<180
%   q2=30*pi/180;         % -45<q2<125
%   q3=30*pi/180;         % -45<q3<180
%   q4=0*pi/180;         %  0<q4<135
%   q5=40*pi/180;         %  -180<q5<180
%   q6=20*pi/180;         %  -90<q6<90
 
%% Matrices de tranformacion por DH
A_BH=[cos(k) 0 -sin(k) 0; 0 1 0 0; sin(k) 0 cos(k) 0; 0 0 0 1];
A_BT=[cos(q1) 0 -sin(q1) 0; sin(q1) 0 cos(q1) 0; 0 -1 0 d1; 0 0 0 1];
A1=A_BH*A_BT;
A2=[cos(q2) 0 sin(q2) a2*cos(q2); sin(q2) 0 -cos(q2) a2*sin(q2); 0 1 0 0; 0 0 0 1];
A3=[cos(q3) 0 -sin(q3) a3*cos(q3); sin(q3) 0 cos(q3) a3*sin(q3); 0 -1 0 d3; 0 0 0 1];
A4=[cos(q4) 0 sin(q4) a4*cos(q4); sin(q4) 0 -cos(q4) a4*sin(q4); 0 1 0 0; 0 0 0 1];
A5=[cos(q5) 0 -sin(q5) 0; sin(q5) 0 cos(q5) 0; 0 -1 0 d5; 0 0 0 1];
AE=[cos(q6) 0 sin(q6) 0; sin(q6) 0 -cos(q6) 0; 0 1 0 0; 0 0 0 1];
BE=[1 0 0 0; 0 1 0 0; 0 0 1 d7; 0 0 0 1];
A6=AE*BE;

%%%Construccion del jacobiano%%%%%%

%%%Para el cálculo del jacobiano se hace uso de la cinemática directa
%%%previamente obtenida. Se definen las siguientes matrices de transformación 
%%%homogénea. 
 
T1=A1;
T2=simplify(A1*A2);
T3=simplify(A1*A2*A3);
T4=simplify(A1*A2*A3*A4);
T5=simplify(A1*A2*A3*A4*A5);
T6=simplify(A1*A2*A3*A4*A5*A6);

%%%De las matrices T_0_i previamente definidas se obtiene ahora los z_i
z0=A_BH(1:3,3);
z1=A1(1:3,3);
z2=T2(1:3,3);
z3=T3(1:3,3);
z4=T4(1:3,3);
z5=T5(1:3,3);

%%%Haciendo también uso de las definiciones de matrices previas se obtiene
%%%los vectores p_i que brindan la posición de i-eslabón.

p0=[0;0;0];
p1=T1(1:3,4);
p2=T2(1:3,4);
p3=T3(1:3,4);
p4=T4(1:3,4);
p5=T5(1:3,4);
pe=T6(1:3,4);
JG=simplify([cross(z0,pe-p0),cross(z1,pe-p1),cross(z2,pe-p2),cross(z3,pe-p3),cross(z4,pe-p4),cross(z5,pe-p5);
             z0              z1              z2              z3              z4              z5           ]);

%% Comprobar el Jacobiano Geom�trico

 %--- Substituci�n de valores 
 %--- Angulos del espacio articular 
 
   q_1=30*pi/180;
   q_2=60*pi/180;
   q_3=90*pi/180;
   q_4=120*pi/180;
   q_5=150*pi/180;
   q_6=180*pi/180;
    
 % Jacobiano num�rico para un punto en el espacio articular (q1,q2,q3,q4,q5,q6)
 JG=subs(JG,[d1, a2, d3, a3, a4, d5, d7, k, q1, q2, q3, q4, q5, q6], [183.73, 58.28, 244.4, 18.76, 18.53, 123.03, 150, 65*pi/180, q_1, q_2, q_3, q_4, q_5, q_6]); 
 JG=double(JG)
