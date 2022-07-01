clc
clear all
close all

% Este script MATLAB permite importar o 
% ficheiro bina'rio produzido pelos programas em C
% e produzir um plot. Pode ser muito util para
% ver se as condicoes fronteira estao corretas,
% para ve se nao ha' erros nas fronteiras dos
% subdominios, etc.

nx=100;
ny=nx;

% alinea = 'a';
% alinea = 'b';
% alinea = 'c';
alinea = 'd';

% fileID = fopen('results_2D.bin');

f = [alinea '/results_' alinea '.bin'];
fileID = fopen(f);

array_MPI = fread(fileID, [ny nx],'double');
fclose(fileID);

L=1;
x=linspace(-L,L,nx);
y=linspace(-L,L,ny);

% Como para as figuras do MATLAB,
% o avan?o numa linha ? um aumento de x 
% (o MATLAB l? em "column-major"),
% para a figura ficar consistente com o 
% programa em C (o nosso programa C escreve 
% em "row-major"), tem que se transpor a matriz.

figure
mesh(x,y,array_MPI')
xlim([-L L])
ylim([-L L])
xlabel('\it{x}')
ylabel('\it{y}')
title('array\_MPI')

i = [alinea,'/img', upper(alinea), '.jpg'];
saveas(gcf,i)


