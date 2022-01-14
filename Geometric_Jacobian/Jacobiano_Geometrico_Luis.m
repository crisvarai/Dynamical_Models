%%%En este scrip se calcula el Jacobiano Geometrico (Producto cruz) del brazo izquierdo del
%%%robot Prometheus, 

syms q1 q2 q3 q4 q5 q6 %Definimos las variables articulares (simbolicos)
syms k d1 a2 a3 d3 a4 d5 d7 %Definimos las variables de parametros fisicos (simbolicos)

%%Matrices de tranformacion por DH
A_BH=[cos(k) 0 -sin(k) 0; 0 1 0 0; sin(k) 0 cos(k) 0; 0 0 0 1];
T1=[cos(q1) 0 -sin(q1) 0; sin(q1) 0 cos(q1) 0; 0 -1 0 d1; 0 0 0 1];
t_New1=A_BH*T1;
T2=[cos(q2) 0 sin(q2) a2*cos(q2); sin(q2) 0 -cos(q2) a2*sin(q2); 0 1 0 0; 0 0 0 1];
T3=[cos(q3) 0 -sin(q3) a3*cos(q3); sin(q3) 0 cos(q3) a3*sin(q3); 0 -1 0 d3; 0 0 0 1];
T4=[cos(q4) 0 sin(q4) a4*cos(q4); sin(q4) 0 -cos(q4) a4*sin(q4); 0 1 0 0; 0 0 0 1];
T5=[cos(q5) 0 -sin(q5) 0; sin(q5) 0 cos(q5) 0; 0 -1 0 d5; 0 0 0 1];
T6=[cos(q6) 0 sin(q6) 0; sin(q6) 0 -cos(q6) 0; 0 1 0 0; 1 0 0 1];
BE=[1 0 0 0; 0 1 0 0; 0 0 1 d7; 0 0 0 1];
TNew7=T6*BE;
T=simplify(t_New1*T2*T3*T4*T5*TNew7); %%%Matriz de transformacion Homogena del marco de referencia global al marco de referencia del End-Effector (Hand)


%%%Construccion del jacobiano%%%%%%

%%%Para el cálculo del jacobiano se hace uso de la cinemática directa
%%%previamente obtenida. Se definen las siguientes matrices de transformación 
%%%homogénea. 
 
T_0_2=t_New1*T2;
T_0_3=t_New1*T2*T3;
T_0_4=t_New1*T2*T3*T4;
T_0_5=t_New1*T2*T3*T4*T5;
T_0_6=T;

%%%De las matrices T_0_i previamente definidas se obtiene ahora los z_i
z0=t_New1(1:3,3);
z1=T_0_2(1:3,3);
z2=T_0_3(1:3,3);
z3=T_0_4(1:3,3);
z4=T_0_5(1:3,3);
z5=T_0_6(1:3,3);

%%%Haciendo también uso de las definiciones de matrices previas se obtiene
%%%los vectores p_i que brindan la posición de i-eslabón.

p0=[0;0;0];
p1=t_New1(1:3,4);
p2=T_0_2(1:3,4);
p3=T_0_3(1:3,4);
p4=T_0_4(1:3,4);
p5=T_0_5(1:3,4);
pe=T_0_6(1:3,4);
JG=simplify([cross(z0,pe-p0),cross(z1,pe-p1),cross(z2,pe-p2),cross(z3,pe-p3),cross(z4,pe-p4),cross(z5,pe-p5);
                          z0,             z1,             z2,             z3,             z4,           z5]);

%% Comprobar el Jacobiano Geom�trico

 %--- Substituci�n de valores 
 %--- Angulos del espacio articular 
 
   q_1=0*pi/180;
   q_2=0*pi/180;
   q_3=0*pi/180;
   q_4=0*pi/180;
   q_5=0*pi/180;
   q_6=0*pi/180;
    
 % Jacobiano num�rico para un punto en el espacio articular (q1,q2,q3,q4,q5,q6)
 JG=subs(JG,[d1, a2, d3, a3, a4, d5, d7, k, q1, q2, q3, q4, q5, q6], [183.73, 58.28, 244.4, 18.76, 18.53, 123.03, 150, 65*pi/180, q_1, q_2, q_3, q_4, q_5, q_6]); 
 JG=double(JG)