%Jonas Nunes
%Modelo teorico de antena ULA
close all;
clc;
clear;
L = 4;              %quantidade de elementos de antena do conjunto              
a = 30*L;           
espacamento = 1- ((a-1)/a);  
x = [-1:espacamento:1];     
theta= acos(x); 

q=[1 3 7 20]  %parametro de diretividade
figure('Name','Feixe do elemento de antena');

for x=1:4
    subplot(2,2,x)
    p_theta = sin(theta).^q(x);    % pattern de radiacao;
    polarplot(theta,abs(p_theta));
    title(['q = ',int2str(q(x))])
end