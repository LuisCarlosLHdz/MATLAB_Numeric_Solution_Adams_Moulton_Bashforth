# MATLAB_Numeric_Solution_Adams_Moulton_Bashforth\
#Example Script -------------------------------------------------\
h = 0.1; %Paso de Integración, Step\
Tlim = 20; %Tiempo Límite, Time limit\
alpha = [0.8,0.8]; %Orden fraccionario, Fractional Order\
n = Tlim/h ; %Nuemro de iteraciones, Number of iterations\
A = 1; %Parametro, Parameter\
B = 3; %Parametro, Parameter\
%----Ecuaciones, Equations, Brusselator System\
syms x z\
eqn1 = A - (B+1)*x + (x^2)*z;\
eqn2 = B*x - (x^2)*z;\
eqns = [eqn1; eqn2]; %Vector de ecuaciones, Vector equations\
%----Condiciones Iniciales, Initial Conditions\
x_t0 = 1.2; \
z_t0 = 2.8;\
con_ini = [x_t0,z_t0]; %Vector de CI, Vecot of IC\
% Adams_Moulton_Bashforth\
%-----Llamar función, Call function "AMB_Method()"\
[y_sol,y_p] = AMB_Method(h,n,alpha,con_ini,eqns);\
%-----Plot soution\
plot(y_sol(1,:),y_sol(2,:))\
----------------------------------------------------------------\
----------------------------------------------------------------\
The function needs the following parameters:\
La función necesita de los siguientes parametros:\
h-> Step\
n-> Number of iterations\
alpha -> Fraciontal order\
con_ini -> Initial conditions\
eqns -> Equations of your system\
----------------------------------------------------------------\
The function return the following vectors:\
La función regresa los siguientes vectores:\
y_sol -> Solution vector\
y_p -> Solution from predictor\
