close all
clear all
clc

%% 
maxit = 5e5;
tolerancia = 1e-6;

nx = 100;
ny = nx;
L = 1;
x = linspace(-L,L,nx);
y = linspace(-L,L,ny);

% h=0.025/L;
% N=2/h+1;
h = 2/nx;
N = nx;
Vold=zeros(N,N);

W = 15/16;

%% 

Vnew = Vold;
for k = 1:maxit
    Vnew = Vold;
    for i = 1:N
        for j = 1:N
            f = 7*sin(2*pi*x(i))*cos(3*pi*x(i))*sin(2*pi*y(j))*cos(3*pi*y(j));

            j1= j + 1;
            j0 = j - 1;
            j2 = j + 2;
            j00 = j - 2;

            if (j0 <= 0)
                j0 = j0+N;
                j00 = j00+N;
            elseif (j1 >= N)
                j1 = j1-N+1;
                j2 = j2-N+1;
            elseif (j00 == 0)
                j00 = N;
            elseif (j2 == N)
                j2 = 1;
            end

            i1= i+1;
            i0 = i-1;
            i2 = i + 2;
            i00 = i - 2;

            if (i0 <= 0)
                i0 = i0+N;
                i00 = i00+N;
            elseif (i1 >= N)
                i1 = i1-N+1;
                i2 = i2-N+1;
            elseif (i00 == 0)
                i00 = N;
            elseif (i2 == N)
                i2 = 1;
            end

            Vnew(i,j)= W/60*(16*Vnew(i0,j)+16*Vnew(i1,j)+16*Vnew(i,j0)+16*Vnew(i,j1)-12*h^2*f + ...
                Vnew(i2,j)+Vnew(i00,j)+Vnew(i,j2)+Vnew(i,j00))+ (1-W)*Vnew(i,j);
        end
    end
    diff = sqrt(sum(sum((Vnew - Vold).^2))) / sqrt(sum(sum(Vnew.^2)));
    if diff < tolerancia
        break;
    end
	Vold=Vnew;
end 

%%
figure;
mesh(x,y,Vnew)
xlim([-L L])
ylim([-L L])
xlabel('\it{x}')
ylabel('\it{y}')
title('array\_MATLAB')

saveas(gcf,"jac_c.jpg")