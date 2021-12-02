%Jonas Nunes
%Gaussian Randomization Procedure
close all;
clc;
clear;
f=3.5e9;          % frequencia
c=3e8;          % velocidade da luz no vacuo

lambda=c/f;     % comprimento de onda
d=lambda/2;     % distancia entre elementos de antenas
etha_0 = 377;   %impedancia do espaco livre;

L = 4 ;                     % qtd de antena no array (linhas)
tam_CodeBook = 8            % tamanho do CodeBook
a = 30*L;                   
espacamento = 1- ((a-1)/a);
x = [-1:espacamento:1];     
theta= acos(x);             

q=3;                        % controla a direcao /ganho do beam formado, quanto maior mais estreito o feixe
p_theta = sin(theta).^q;    % pattern de radiacao do elemento de antena simples
e_phi = 0;                  % considerado nulo pois e o eixo de formacao do conjunto ULA
                           
e_theta = sqrt(p_theta) .* exp(j * (2*pi*cos(theta))/d .* [0:L-1]');

for pattern=1:1:size(theta,2)
vetor_M(:,:,pattern) = e_theta(:,pattern) * ctranspose(e_theta(:,pattern)) + e_phi * ctranspose(e_phi);
end

M = vetor_M(:,:,pattern);
pattern = 1;
k=1;
w0 = randn(L,k) + j*randn(L,k);
%[w0,~] = eigs(vetor_M(:,:,pattern),1);
W0 = w0 * ctranspose(w0);
%autovalor = lambda
[autovetores,autovalores] = eigs(W0);

tam_CodeBook=4;
%notacao padrao i = linha, k = coluna da matriz M (nesse exemplo 4x4)
%aqui ira gerar um M para cada pattern que e de acordo com a variacao do
%angulo theta, logo tera uma matriz de resposta do campo eletrico para cada
%um dos patterns possiveis.
for pattern=1:1:size(theta,2)
vetor_M(:,:,pattern) = e_theta(:,pattern) * ctranspose(e_theta(:,pattern)) + e_phi * ctranspose(e_phi);
end

pattern=1; %como e formacao de beam unico, selecionei qualquer pattern, poderia colocar um for para fazer cada uma das possibilidades
K=tam_CodeBook;
W = zeros(L,K);

M = vetor_M(:,:,pattern);
%notacao padrao i = linha, k = coluna da matriz M (nesse exemplo 4x4)
w=[1;j;-1;1]; %Peso qualquer aleatorio para um feixe  qualquer, ou seja
% so possui uma coluna, ja que por cima dele que os outros patterns serao
% formados (outras colunas)

    somatorio=0;
    parcial = 0;
    G_antigo = 0;
    G = 1;
    vetor_g=0;
    vetor_w = zeros(4,1);
    i=1;
    %vai gerando os pesos aleatorios e comparando eles com o proximo
    % como esta no algoritimo 2
    while G>G_antigo
        for ii = i:K%size(theta,2)
            i= mod(ii, L)+1;
        for k = 1:size(M,2)
            if k==i
                parcial = 0;
            else
                parcial = M(i,k) * w(i);
            end
                somatorio = somatorio + parcial;
        end
        w(i) = 1/sqrt(L) * exp(j*angle(somatorio)); %equacao 17

      G = 2*pi/etha_0 * ctranspose(w)*M *w;        %equacao 7
      vetor_g(length(vetor_g)+1) = G;
      vetor_w(:,size(vetor_w,2)+1) = w;
      G_antigo = G;
        end
    end
    
    %remover a primeira coluna pois ela e composta apenas de zeros
    vetor_g(:,1) = [];
    vetor_w(:,1) = [];