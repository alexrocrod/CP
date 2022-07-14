clc
clear all
close all

%% 

nx=100;
ny=nx;

% alinea = 'a';
% alinea = 'b';
% alinea = 'c';
alinea = 'd';


f = ['../results_' alinea '.bin'];

fileID = fopen(f);

array_MPI = fread(fileID, [ny nx],'double');
fclose(fileID);

fexact = ['vnewMat_' alinea '.mat'];
load(fexact)

L=1;
x=linspace(-L,L,nx);
y=linspace(-L,L,ny);

figure
mesh(x,y,array_MPI')
xlim([-L L])
ylim([-L L])
xlabel('\it{x}')
ylabel('\it{y}')
title('array\_MPI')
fprintf("MSE: %d\n",getMSE(Vnew,array_MPI',nx));

i = ['../img', upper(alinea), '.jpg'];

saveas(gcf,i)
%%

figure
mesh(x,y,Vnew)
xlim([-L L])
ylim([-L L])
xlabel('\it{x}')
ylabel('\it{y}')
title('array\_Matlab')


%%

function MSE = getMSE(Vnew,array_MPI,nx)

    N2 = nx*nx;
    
    MSE = 1/N2 * sum((array_MPI-Vnew).^2,'all');
end




