% Jonas Nunes
% Visualizador de CodeBook
clear all;
close all;

f=3.5e9; ; %frequencia
c=3e8; % velocidade da luz no vacuo
lambda = c/f; %comprimento de onda
d = lambda/2;  % distancia dos elementos de antena
etha_0 = 377;   %impedancia do espaco livre;

% Numero de elementos de antenas
M = 4;
% Numero de patterns (codewords)
K =4;

W = [1,1,1,1; ...
    -1,-j,1,j; ...
    1,-1,1,-1;...
    -1,j,1,-j];


theta = 0:0.01:2*pi; %resolucao da formacao do diagrama de radiacao, ``quantidade de direcoes do CodeBook"


K_theta = length(theta);
arrayFactor = zeros(K,K_theta);
for k_theta=1:K_theta
    for k=0:K-1
        for m=0:M-1
            arrayFactor(k+1,k_theta) = arrayFactor(k+1,k_theta) + W(m+1,k+1)* exp(i*2*pi*m*(d/lambda)*cos(theta(k_theta)));
        end
    end
end

%juntar os maximos de cada para criar um beam com o maximo do conjunto em
%cada direcao
for i = 1:size(arrayFactor,2)
arrayFactorTotal(i) = max(arrayFactor(:,i));
end

% Plot de cada um dos patterns do CodeBook
figure;
stringPlot={'b','r','g','c','m','y','k','w'};
for k=1:K
    polar(theta,abs(arrayFactor(k,:)),stringPlot{k})
    if k==1 
        hold on;
    end
    stringLegend{k}= ['id-pattern-' num2str(k)];
end

hold off
legend(stringLegend)

figure,polar(theta,abs(arrayFactorTotal),stringPlot{k})%plota o pattern conjunto
legend(``ArrayFactor Total'')
xlabel('\theta (rad)')

figure, plot(theta-pi,abs(arrayFactor)') %plot de cada beam 2D
legend(stringLegend)
xlabel('\theta (rad)')
