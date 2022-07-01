close all
clear all
clc

L = 1;
h=0.025/L;
N=2/h+1;
Vold=zeros(N,N);

Ncent = round(N/4);
Vold(Ncent+1:end-Ncent,Ncent+1:end-Ncent) = 1;
condicao = logical(Vold);
tolerancia = 1e-7;

N_max = 10000;
Vnew = Vold;
for k = 1:N_max
    Vnew = Vold;
    for i = 2:N-1
        for j = 2:N-1
            f(i,j) = -2*(2-x(i)^2-y(j)^2);           
            Vnew(i,j)= 0.25*(Vnew(i+1,j)+Vnew(i-1,j)+Vnew(i,j+1)+Vnew(i,j-1)-h^2*f(i,j));
        end
    end
    diff = sqrt(sum(sum((Vnew - Vold).^2))) / sqrt(sum(sum(Vnew.^2)));
    if diff < tolerancia
        break;
    end
	Vold=Vnew;
end 