%Jonas Nunes
close all;
clc;
clear;

f=3.5e9;          % frequencia
c=3e8;          % velocidade da luz no vacuo
K=16;            % tamanho do CodeBook (colunas)
lambda=c/f;     % comprimento de onda
d=lambda/2;     % distancia entre elementos de antenas
etha_0 = 377;   %impedancia do espaco livre;

L = 4 ;                     % qtd de antena no array (linhas)
a = 30*L;                   
espacamento = 1- ((a-1)/a);  
x = [-1:espacamento:1];     % (x=-1,...,1)
theta= acos(x);             
q=3;                        % controla a direcao /ganho do beam formado, quanto maior mais estreito o feixe
p_theta = sin(theta).^q;    % pattern de radiacao do elemento de antena simples
e_phi = 0;                  % considerado nulo pois e o eixo de formacao do conjunto ULA
                            %e_theta descrito na equacao 33
e_theta = sqrt(p_theta) .* exp(j * (2*pi*cos(theta))/d .* [0:L-1]');

figure('Name','Feixe do elemento de antena');
polarplot(theta,abs(sin(theta).^q));

%vetor_M apresenta todas as matrizes M (LxK) x Angulos de theta (pattern)
for angulo=1:1:size(theta,2)
vetor_M(:,:,angulo) = e_theta(:,angulo) * ctranspose(e_theta(:,angulo)) + e_phi * ctranspose(e_phi);
end

%geracao dos beams individuais
for angulo=1:1:size(theta,2)
%notacao padrao i = linha, k = coluna da matriz M (nesse exemplo 4x4)
w2 = randn(L,1) + j*randn(L,1);
% so possui uma coluna, ja que por cima dele que os outros patterns serao
% formados (outras colunas)

M = vetor_M(:,:,angulo)
    somatorio=0;
    parcial = 0;
    G_antigo = 0;
    G = 1;
    vetor_g=0;
    vetor_w = zeros(L,1);
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
                parcial = M(i,k) * w2(i);
            end
                somatorio = somatorio + parcial;
        end
        w2(i) = 1/sqrt(L) * exp(j*angle(somatorio)); %equacao 17

      G = 2*pi/etha_0 * ctranspose(w2)*M *w2;        %equacao 7
      vetor_g(length(vetor_g)+1) = G;
      vetor_w(:,size(vetor_w,2)+1) = w2;
      G_antigo = G;
        end
    end
    
    %remover a primeira coluna pois ela e composta apenas de zeros
    vetor_g(:,1) = [];
    vetor_w(:,1) = [];
    %vetor_wzao apresenta todos esses pesos calculados (LxK) x Angulos de theta (pattern)
    vetor_wzao(:,:,angulo)= vetor_w(:,:);
end

for pattern = 1:size(theta,2) %pattern 8 pra 4 legal
    %para saber o que e cada angulo em grau e so usar rad2deg(theta(pattern))
    M = vetor_M(:,:,pattern)
[autovetores(:,i),autovalormax(:,i)] = eigs(M,1);
% w e o autovetor correspondente ao maior autovalor

w(:,i) = autovetores(:,i);
S= 2*pi/etha_0 * max(ctranspose(w(:,i))*M*w(:,i));
vetor_S(i) = S;
i=i+1;
end

%pegar os indices dos pesos que geram as maiores ganho maximo de cada
%CodeWord (S = ganho maximo de um CodeWord na esfera teorica)
[~,indices] = sort(vetor_S,'descend')

clear W;
for i=1:size(indices,2)
W(:,i) = w(:,indices(i));
end

K_theta = length(theta);
arrayFactor = zeros(K,K_theta);
theta_max = zeros(size(arrayFactor,1),2);

for k_theta=1:K_theta
    for k=0:K-1
        for l=0:L-1
            arrayFactor(k+1,k_theta) = arrayFactor(k+1,k_theta) + W(l+1,k+1)* exp(j*2*pi*l*(d/lambda)*cos(theta(k_theta)));
        end
    end
end

theta_max = zeros(size(arrayFactor,1),1);
%arrayFactor = arrayFactor/sqrt(2); % pra ficar igual aos outros com maximo =1.5

%juntar todos os beam em um so
for i = 1:size(arrayFactor,2)
arrayFactorTotal(i) = max(arrayFactor(:,i));
end

%Normalizacao
arrayFactorTotal = arrayFactorTotal/max(arrayFactorTotal)
for k=1:K
    arrayFactor(k,:) = arrayFactor(k,:)/max(arrayFactor(k,:))
end

% Plot graphs
figure;
stringPlot={'b','r','g','c','m','y','k','b','r','g','c','m','y','k'};
for k=1:K
    polar(theta,abs(arrayFactor(k,:)),stringPlot{k})
    hold on;
    %para plotar a parte inferior do beam que e espelhada, uma vez que o
    %modelo so calculou de 0 ate 180.
    title(['CodeBook - Greedy K =',num2str(k),' L =',num2str(L)]);
    if k==1 
        hold on;
    end
    stringLegend{k}= ['id-pattern-' num2str(k)];
end

%plotar a parte inferior do diagrama de radiacao.
for k=1:K
    hold on;
    polar(-theta,abs(arrayFactor(k,:)),stringPlot{k})
end
hold off
legend(stringLegend)

figure
hold on
for k=1:K
plot(rad2deg(theta-pi),abs(arrayFactor(k,:))',stringPlot{k})
title(['CodeBook - Greedy K =',num2str(k),' L =',num2str(L)]);
xlabel('\theta (graus)')
end

%perda maxima = diferenca entre valor maximo e minimo
max_loss = max(abs(arrayFactorTotal)) - min(abs(arrayFactorTotal))

legend(stringLegend)
hold off
%theta_max(:,2) = radtodeg(theta_max(:,1))
%ver de remover
figure;
polar(theta,abs(arrayFactorTotal),stringPlot{1})
%para plotar a parte inferior do beam que e espelhada, uma vez que o
    %modelo so calculou de 0 ate 180.
hold on;
polar(-theta,abs(arrayFactorTotal),stringPlot{1})
title(['CodeBook - Greedy K =',num2str(k),' L =',num2str(L)]);
xlabel('\theta (graus)')
annotation('textbox',[.80 .7 .2 .2],'String',['Perda Max: ' num2str(max_loss)],'EdgeColor','none')

clear U20;
clear U50;
clear pattern_id;
clear Umean;
clear theta_max
%Tabela das informacoes
for i=1:size(arrayFactor,1)
[maximo,idx] = max(abs(arrayFactor(i,:)));
theta_max(i) = rad2deg(theta(idx));
Umean(i) = mean(db(arrayFactor(i,:)));
U20(i) = prctile(db(arrayFactor(i,:)),20);
U50(i) = prctile(db(arrayFactor(i,:)),50);
pattern_id(i) = strcat(``Pattern - `` , int2str(i));
end
pattern_id(size(pattern_id,2)+1) = ``CodeBook''
theta_max(size(theta_max,2)+1) = ``N/A"
UmeanTotal = mean(db(arrayFactorTotal));
U20Total = prctile(db(arrayFactorTotal),20);
U50Total = prctile(db(arrayFactorTotal),50);
Umean(size(Umean)+1) = UmeanTotal
U20(size(U20)+1) = U20Total
U50(size(U50)+1) = U50Total
theta_max = theta_max'
Umean = Umean'
U20 = U20'
U50 = U50'
pattern_id = pattern_id'
table(pattern_id,theta_max,Umean,U50,U20)