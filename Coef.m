function [a, b] = Coef(h,n,alpha)
%Coeficientes del sistema numérico Adams Moulton V2.0
%M.Cs. Luis Carlos Lujano Hernández
%   Función Coeficeintes a, b son vectores de tamaño n+1
[~,Num_alpha] = size(alpha);
b = zeros(Num_alpha,n+1);
a = zeros(Num_alpha,n+1);
for i=1:Num_alpha
    for j=1:n+1
        %----a_j,n+1
        if j==1
            a(i,j) = ((n)^(alpha(1,i) + 1)) - ( (n - alpha(1,i))*((n + 1)^(alpha(1,i))) );
        elseif j>=2 && j<=(n+1)
            a(i,j) = ((n - (j-1) + 2)^(alpha(1,i) + 1)) + ((n - (j-1))^(alpha(1,i) + 1)) - (2*((n - (j-1) + 1)^(alpha(1,i) + 1)));
        else
            a(i,j) = 1;
        end
        %----b_j,n+1
        b(i,j) = (((h)^(alpha(1,i)))/alpha(1,i))*( ( (n + 1 - (j-1))^(alpha(1,i)) ) - ( (n - (j-1))^(alpha(1,i)) ) );
    end
end
end