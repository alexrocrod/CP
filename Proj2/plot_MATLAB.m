clc
clear all
close all

nx=100;
ny=nx;

% alinea = 'a';
% alinea = 'b';
% alinea = 'c';
alinea = 'd';

C = false;
% C = true;

f = [alinea '/results_' alinea '.bin'];

if (C)
    f = ['d/results_' alinea '_C.bin'];
end

fileID = fopen(f);

array_MPI = fread(fileID, [ny nx],'double');
fclose(fileID);

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

i = [alinea,'/img', upper(alinea), '.jpg'];
if (C)
    i = [alinea,'/img', upper(alinea), '_C.jpg'];
end
saveas(gcf,i)


