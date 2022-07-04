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


%% 

Vnew = Vold;
for k = 1:maxit
%     Vnew = Vold;
    for i = 1:N
        for j = 1:N
            f = 7*sin(2*pi*x(i))*cos(3*pi*x(i))*sin(2*pi*y(j))*cos(3*pi*y(j));

            i1= i+1;
            i0 = i-1;
            j1= j+1;
            j0 = j-1;

            if (j0 == 0)
                j0 = N;
            end
            if (j1 == N+1)
                j1 = 1;
            end

            if (i0 == 0)
                i0 = N;
            end
            if (i1 == N+1)
                i1 = 1;
            end
                 
            Vnew(i,j)= 0.25*(Vnew(i1,j)+Vnew(i0,j)+Vnew(i,j1)+Vnew(i,j0)-h^2*f);
        end
    end
    diff = sqrt(sum(sum((Vnew - Vold).^2))) / sqrt(sum(sum(Vnew.^2)));
    if diff < tolerancia
        break;
    end
	Vold = Vnew;
end 

%%
figure;
mesh(x,y,Vnew)
xlim([-L L])
ylim([-L L])
xlabel('\it{x}')
ylabel('\it{y}')
title('array\_MATLAB')

saveas(gcf,"jac_d.jpg")