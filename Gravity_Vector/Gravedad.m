% compensación de gravedad Prometheus 

clc
close all
clear all

%----- Variables 
% parametros del robot
syms k d1 a2 a3 d3 a4 d5 d7
% variables articulares 
syms q1 q2 q3 q4 q5 q6 
% masas  
syms m1 m2 m3 m4 m5 m6
%-- Aceleración de la gravedad 
syms g
%-- Cordenada del cp desde la referencia Oi, i=1,2,3,4,5,6 
syms x1 y1 z1
syms x2 y2 z2
syms x3 y3 z3
syms x4 y4 z4
syms x5 y5 z5
syms x6 y6 z6

%----- BRAZO IZQUIERDO 
%----- Matrices Homogeneas
TLBH=[cos(k) 0 -sin(k) 0; 0 1 0 0; sin(k) 0 cos(k) 0; 0 0 0 1];
TL01=[cos(q1) 0 -sin(q1) 0; sin(q1) 0 cos(q1) 0; 0 -1 0 d1; 0 0 0 1];
TL12=[cos(q2) 0 sin(q2) a2*cos(q2); sin(q2) 0 -cos(q2) a2*sin(q2); 0 1 0 0; 0 0 0 1];
TL23=[cos(q3) 0 -sin(q3) a3*cos(q3); sin(q3) 0 cos(q3) a3*sin(q3); 0 -1 0 d3; 0 0 0 1];
TL34=[cos(q4) 0 sin(q4) a4*cos(q4); sin(q4) 0 -cos(q4) a4*sin(q4); 0 1 0 0; 0 0 0 1];
TL45=[cos(q5) 0 -sin(q5) 0; sin(q5) 0 cos(q5) 0; 0 -1 0 d5; 0 0 0 1];
TL56=[cos(q6) 0 sin(q6) 0; sin(q1) 0 -cos(q6) 0; 0 1 0 0; 0 0 0 1];
TL67=[1 0 0 0; 0 1 0 0; 0 0 1 d7; 0 0 0 1];

%----- Matriz de orientación-posición en el eje de rotación de cada actuador 
A1=simplify(TLBH*TL01);
A2=simplify(TLBH*TL01*TL12);
A3=simplify(TLBH*TL01*TL12*TL23);
A4=simplify(TLBH*TL01*TL12*TL23*TL34);
A5=simplify(TLBH*TL01*TL12*TL23*TL34*TL45);
A6=simplify(TLBH*TL01*TL12*TL23*TL34*TL45*TL56*TL67);

%---- Posiciónes 
O10=[A1(1,4); A1(2,4); A1(3,4)];
O20=[A2(1,4); A2(2,4); A2(3,4)];
O30=[A3(1,4); A3(2,4); A3(3,4)];
O40=[A4(1,4); A4(2,4); A4(3,4)];
O50=[A5(1,4); A5(2,4); A5(3,4)];
O60=[A6(1,4); A6(2,4); A6(3,4)];

%---- Rotaciones
R10=[A1(1,1) A1(1,2) A1(1,3); A1(2,1) A1(2,2) A1(2,3); A1(3,1) A1(3,2) A1(3,3)];
R20=[A2(1,1) A2(1,2) A2(1,3); A2(2,1) A2(2,2) A2(2,3); A2(3,1) A2(3,2) A2(3,3)];
R30=[A3(1,1) A3(1,2) A3(1,3); A3(2,1) A3(2,2) A3(2,3); A3(3,1) A3(3,2) A3(3,3)];
R40=[A4(1,1) A4(1,2) A4(1,3); A4(2,1) A4(2,2) A4(2,3); A4(3,1) A4(3,2) A4(3,3)];
R50=[A5(1,1) A5(1,2) A5(1,3); A5(2,1) A5(2,2) A5(2,3); A5(3,1) A5(3,2) A5(3,3)];
R60=[A6(1,1) A6(1,2) A6(1,3); A6(2,1) A6(2,2) A6(2,3); A6(3,1) A6(3,2) A6(3,3)];

%---- Posición del centro de masa respecto al sistema de referencia inercial 
%---- Dados por Alex 

   cm1=[-106.41; -0.92; 49.68];
   cm2=[-199.44; 1.53; 151.77];
   cm3=[-294.44; 41.45; 222.04];
   cm4=[-399.42; -0.01; 291.26];
   cm5=[-458.11; -0.39; 318.81];
   cm6=[-600.99; -6.07; 371.74];

%---- posición p_i del Cm respecto al sistema de referencia i
  P1_s=-transpose(R10)*O10 + transpose(R10)*cm1;
  P2_s=-transpose(R20)*O20 + transpose(R20)*cm2;
  P3_s=-transpose(R30)*O30 + transpose(R30)*cm3;
  P4_s=-transpose(R40)*O40 + transpose(R40)*cm4;
  P5_s=-transpose(R50)*O50 + transpose(R50)*cm5;
  P6_s=-transpose(R60)*O60 + transpose(R60)*cm6;

P1=[x1; y1; z1];
P2=[x2; y2; z2];
P3=[x3; y3; z3];
P4=[x4; y4; z4];
P5=[x5; y5; z5];
P6=[x6; y6; z6];

% Posición del centro de masa en función de las variables articulares qi, i=1,2,3,4,5,6 

Pcm1=O10 + R10*P1;
Pcm2=O20 + R20*P2;
Pcm3=O30 + R30*P3;
Pcm4=O40 + R40*P4;
Pcm5=O50 + R50*P5;
Pcm6=O60 + R60*P6;

%---- Energía potencial, la gravedad unicamente infiere en la componente z 
U1=m1*g*Pcm1(3);
U2=m2*g*Pcm2(3);
U3=m3*g*Pcm3(3);
U4=m4*g*Pcm4(3);
U5=m5*g*Pcm5(3);
U6=m6*g*Pcm6(3);

%---- Energía potencial total 
 U=U1+U2+U3+U4+U5+U6;
 
%--- Gradiente 
 
 DU_q1=diff(U,q1);
 DU_q2=diff(U,q2);
 DU_q3=diff(U,q3);
 DU_q4=diff(U,q4);
 DU_q5=diff(U,q5);
 DU_q6=diff(U,q6);
 
 %--- Vector de gravedad 
 G=[DU_q1; DU_q2; DU_q3; DU_q4; DU_q5; DU_q6]
 
 %G=simplify([DU_q1; DU_q2; DU_q3; DU_q4; DU_q5; DU_q6]); % simplicación "tarda mucho en compilar"


%--- substitución de valores 

   %----- Parámetros del brazo, distancias [m]
   d_1=183.73/1000;
   a_2=58.28/1000;
   d_3=244.4/1000; 
   a_3=(13.76+5)/1000;
   a_4=18.53/1000;
   d_5=(78.53+44.5)/1000;
   d_7=150/1000;

   %--- Rotación inicial del hombro
   k_0=((65*pi)/180);

   %----- masas [kg]
   m_1=1929.76/1000;
   m_2=992.51/1000;
   m_3=1126.28/1000;
   m_4=540.00/1000;  
   m_5=566.31/1000;
   m_6=469.22/1000;

   %---- Posiciones P_i del centro de masa respecto al sistema i [m]

    P1=[0.0546; 66.2941; -0.9200]/1000;
    P2=[-5.0167; 1.5300; 61.1648]/1000;
    P3=[ -0.2391; 67.4386; 41.4500]/1000;
    P4=[-0.4010; -0.0100; 56.9593]/1000;
    P5=[-0.2357; 1.2364; -39.0000]/1000;
    P6=[-12.6485; -6.0700; 0.6260]/1000;

    %---- Aceleración de la gravedad [m/s^2]
    g_0=9.81;
        
   %--- Angulos del espacio articular 
   q_1=0*pi/180;         % -20<q1<180
   q_2=0*pi/180;         % -45<q2<125
   q_3=0*pi/180;         % -45<q3<180
   q_4=0*pi/180;         %  0<q4<135
   q_5=0*pi/180;         %  -180<q5<180
   q_6=0*pi/180; 


%---- Vector de gravedfad para un punto en el espacio articular (q1,q2,q3,q4,q5,q6)
G_num=subs(G,[d1, a2, d3, a3, a4, d5, d7, k, q1, q2, q3, q4, q5, q6, x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, z5, x6, y6, z6, m1, m2, m3, m4, m5 m6, g], [d_1, a_2, d_3, a_3, a_4, d_5, d_7, k_0, q_1, q_2, q_3, q_4, q_5, q_6, P1(1), P1(2), P1(3), P2(1), P2(2), P2(3), P3(1), P3(2), P3(3), P4(1), P4(2), P4(3), P5(1), P5(2), P5(3), P6(1), P6(2), P6(3), m_1, m_2, m_3, m_4, m_5 m_6, g_0]); 
G_num=double(G_num)

P1_num=subs(P1_s,[d1,k,q1],[d_1,k_0,q_1]);
P1_num=double(P1_num)

P2_num=subs(P2_s,[d1,k,q1,q2,a2],[d_1,k_0,q_1,q_2,a_2]);
P2_num=double(P2_num)

P3_num=subs(P3_s,[d1,k,q1,q2,a2,d3,a3,q3],[d_1,k_0,q_1,q_2,a_2,d_3,a_3,q_3]);
P3_num=double(P3_num)

P4_num=subs(P4_s,[d1,k,q1,q2,a2,d3,a3,q3,a4,q4],[d_1,k_0,q_1,q_2,a_2,d_3,a_3,q_3,a_4,q_4]);
P4_num=double(P4_num)
  
   
