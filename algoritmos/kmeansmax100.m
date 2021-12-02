%Jonas Nunes
%Algoritiimo Kmeans maior dentre 100 iteracoes
close all;
clc;
clear;

f=3.5e9;          % frequencia
c=3e8;          % velocidade da luz no vacuo
K=5;            % tamanho do codebook (colunas)
tam_codebook=8
lambda=c/f;     % comprimento de onda
d=lambda/2;     % distancia entre elementos de antenas
etha_0 = 377;   %impedancia do espaco livre;

L = 4 ;                     % qtd de antena no array (linhas)
a = 30*L;                   % 
espacamento = 1- ((a-1)/a);  % 
x = [-1:espacamento:1];     % 
theta= acos(x);             % 
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
M = vetor_M(:,:,angulo);

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

i=1;
for pattern = 1:size(theta,2) %pattern 8 pra 4 legal
[autovetores(:,i),autovalormax(:,i)] = eigs(vetor_M(:,:,pattern),1);
i=i+1;
end

% os pesos otimos serao os autovetores atrelados aos maximos autovalores
W = autovetores;
K_theta = length(theta);

arrayFactor = zeros(K_theta,K_theta);
theta_max = zeros(size(arrayFactor,1),2);

 
for k_theta=1:K_theta
    for k=0:K_theta-1
        for l=0:L-1
            arrayFactor(k+1,k_theta) = arrayFactor(k+1,k_theta) + W(l+1,k+1)* exp(j*2*pi*l*(d/lambda)*cos(theta(k_theta)));
        end
    end
end

for jj = 1:100 %%quantidade de iteracoes do kmeans
theta_max = zeros(size(arrayFactor,1),1);

controle=0
Umean_antigo=0

while controle == 0
for k=1:size(arrayFactor,1)
[ganhomax(k,1),index(k,1)] = max(arrayFactor(k,:));
angulo(k,1) = theta(index(k,1));
angulo(k,2) = rad2deg(angulo(k,1));
ganhomax(k,2) = db(ganhomax(k,1));
end

coordenadas = [angulo(:,2),ganhomax(:,2)];
[idx,C] = kmeans(coordenadas,K); %idx tem a cada centorid q pertence e C tem os centroids

for i =1:size(C,1)
distancias(:,i) = abs(coordenadas(:,1) - C(i,1));
[~,indexangulominimo(i)] = min(distancias(:,i));
end

i=1;
 for k=indexangulominimo
	arrayFactorKmeans(jj,i,:)= arrayFactor(k,:);
	 i=i+1;
 end
 
for i = 1:size(arrayFactorKmeans,3)
arrayFactorKmeansTotal(jj,i) = max(arrayFactorKmeans(jj,:,i));
end
 
 Umean = mean(db(arrayFactorKmeansTotal(jj,:)));

 if Umean_antigo < Umean
     Umean_antigo = Umean;
 else
     controle = 1;
     Umean = Umean_antigo;
 end
 
end 
 
 i=1;   
Umean_otimo(jj) = Umean;
 txt = ['Umean: ' num2str(Umean)];
end

 figure;
stringPlot={'b','r','g','c','m','y','k','b','r','g','c','m','y', 'k','b','r','g','c','m','y','k','b','r','g','c','m','y', 'k','b','r','g','c','m','y','k','b','r','g','c','m','y','k'};
hold on;
for i=1:size(C,1)
scatter(C(i,1),C(i,2),[stringPlot{i},'*']);
end

for i = 1: size(C,1)
stringLegend{i}= ['Cluster - ' num2str(i)];
end

[~,indice] = max(Umean_otimo);  

arrayFactorKmeansmax = arrayFactorKmeans(indice,:,:);
arrayFactorKmeansmax = reshape(arrayFactorKmeansmax,[size(arrayFactorKmeansmax,2), size(arrayFactorKmeansmax,3)]);
arrayFactorKmeansTotalmax = arrayFactorKmeansTotal(indice,:,:);

for i=1:size(coordenadas,1)
    title('divisao clusters')
scatter(coordenadas(i,1),coordenadas(i,2),stringPlot{idx(i)});
ylabel('Ganho Maximo (dB)');
xlabel('Angulo de Maximo Ganho (\theta)')
end
hold off;
legend(stringLegend)

clear indexangulominimo;
clear distancias;
%proximidade para estimar o kmeans, estipular o centroids
 
%normalizar
arrayFactorKmeansmax = arrayFactorKmeansmax/max(max(arrayFactorKmeansmax));
arrayFactorKmeansTotalmax = arrayFactorKmeansTotalmax/max(max(arrayFactorKmeansTotalmax));

figure;
 for i=1:size(arrayFactorKmeansmax,1)
     polarplot(theta,abs(arrayFactorKmeansmax(i,:)),stringPlot{i});
     title({'CodeBook - Kmeans'});
     hold on;
     stringLegend{i}= ['id-pattern-' num2str(i)];
 end

%plotar parte inferior
 for i=1:size(arrayFactorKmeansmax,1)
polarplot(-theta,abs(arrayFactorKmeansmax(i,:)),stringPlot{i});
 end
 legend(stringLegend)
 
%plot retangular
 figure
hold on
for k=1:K
plot(rad2deg(theta-pi),abs(arrayFactor(k,:))',stringPlot{k})
title(['CodeBook - Kmeans K =',num2str(k),' L =',num2str(L)]);
xlabel('\theta (graus)')
stringLegend{i}= ['id-pattern-' num2str(i)];
end
legend(stringLegend)
 
max_loss = max(abs(arrayFactorKmeansTotalmax)) - min(abs(arrayFactorKmeansTotalmax)) 
figure;
polarplot(theta,abs(arrayFactorKmeansTotalmax),'r');
title({'CodeBook - Kmeans'});
annotation('textbox',[.80 .7 .2 .2],'String',['Perda Max: ' num2str(max_loss)],'EdgeColor','none')

%plotar parte inferior
hold on;
polarplot(-theta,abs(arrayFactorKmeansTotalmax),'r');
hold off;

%apagar algumas variaveis ja utilizadas anteriormente
clear pattern_id;
clear Umean;
clear theta_max
%Tabela das informacoes
for i=1:size(arrayFactorKmeansmax,1)
[maximo,idx] = max(abs(arrayFactorKmeansmax(i,:)));
theta_max(i) = rad2deg(theta(idx));
Umean(i) = mean(db(arrayFactorKmeansmax(i,:)));
U20(i) = prctile(db(arrayFactorKmeansmax(i,:)),20);
U50(i) = prctile(db(arrayFactorKmeansmax(i,:)),50);
pattern_id(i) = strcat("Pattern - " , int2str(i));
end
pattern_id(size(pattern_id,2)+1) = "CodeBook"
theta_max(size(theta_max,2)+1) = "N/A"
UmeanTotal = mean(db(arrayFactorKmeansTotalmax));
U20Total = prctile(db(arrayFactorKmeansTotalmax),20);
U50Total = prctile(db(arrayFactorKmeansTotalmax),50);
Umean(size(Umean)+1) = UmeanTotal
U20(size(U20)+1) = U20Total
U50(size(U50)+1) = U50Total
theta_max = theta_max'
Umean = Umean'
U20 = U20'
U50 = U50'
pattern_id = pattern_id'
table(pattern_id,theta_max,Umean,U50,U20)
