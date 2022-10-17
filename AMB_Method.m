function [y_sol,y_p] = AMB_Method(h,n,alpha,con_ini,eqns)
%----------------------------------------------------------------------
%----------------------------------------------------------------------
%-------Metodo predictivo-correctivo: Adams_Moulton_Bashforth----------
%-------M.Cs Luis Carlos Lujano Hernández------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------
%-------h       -> paso de integración---------------------------------
%-------n       -> número de pasos-------------------------------------
%-------alpha   -> orden fraccionario [alpha_x,alpha_y,alpha_z]--------
%-------con_ini -> Vector de condiciones iniciales [x0,y0,z0]----------
%-------eqns    -> Vector de ecuaciones [eqn1;eqn2;eqn3]---------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------
[~,Num_var] = size(con_ini); %Número de variables
[Num_eqns,~] = size(eqns); %Número de variables
[a,b] = Coef(h,n,alpha); %Coeficinetes a, b
% -------
var = symvar(eqns); %Variables en la ecuación
ecu_eval = subs(eqns,var,con_ini); %Evalua las ecuaciones en condiciones iniciales
%--------------------------------------------
sum_p  = zeros(Num_eqns,1); %Memoria de suma para el predictivo
y_p    = zeros(1,n+1); %Memoria de la solución del predictivo
sum_cp = zeros(Num_var,1); %Memoria de suma para correctivo-predictivo
y_sol  = zeros(Num_eqns,n+1); %Vector de solución
for m=1:n+1
    if m==1
        %--------Solucion de condicion incial
        for i=1:Num_eqns
            sum_p(i,1)  = (b(i,m)*ecu_eval(i,1)); %Resolver la sumatoria de Adams-Bashforth
            y_p(i,m)    = con_ini(1,i) + ((1/gamma(alpha(1,i)))*sum_p(i,1)); %Obtener la solución Adams-Bashforth
            sum_cp(i,1) = (a(i,m)*ecu_eval(i,1)); %Resolver la sumatoria de Adams-Moulton
            y_sol(i,m)  = con_ini(1,i) + (((h^alpha(1,i))/(gamma(alpha(1,i) + 1)))*y_p(i,m)) + ( ((h^alpha(1,i))/(gamma(alpha(1,i) + 2)))*sum_cp(i,1) ); %Obtener la solución Adams-Moulton
        end
    else
        %-------Solucion para las siguientes
        Nuevo_val = y_sol(:,m-1)'; %Asignar Adams-Moulton anterior a las variables
        ecu_eval = subs(eqns,var,Nuevo_val); %Evalua las ecuaciones en Adams-Moulton
        for i=1:Num_eqns
            sum_p(i,1) = (b(i,m)*ecu_eval(i,1)) + sum_p(i,1); %Resolver la sumatoria de Adams-Bashforth más las soluciones anteriores
            y_p(i,m) = con_ini(1,i) + ((1/gamma(alpha(1,i)))*sum_p(i,1)); %Obtener la solución Adams-Bashforth
            sum_cp(i,1) = (a(i,m)*ecu_eval(i,1)) + sum_cp(i,1); %Resolver la sumatoria de Adams-Moulton más las soluciones anteriores
            y_sol(i,m) = con_ini(1,i) + (((h^alpha(1,i))/(gamma(alpha(1,i) + 1)))*y_p(i,m)) + ((h^alpha(1,i))/(gamma(alpha(1,i) + 2)))*sum_cp(i,1); %Obtener la solución Adams-Moulton
        end
        %---------------
    end
end
end