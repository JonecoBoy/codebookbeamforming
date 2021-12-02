clc;
close all;
clear all;
M = 2; %numero de antenas transmissoras
N = M; %numero de antenas receptoras
ITER = 1000; %numero de tentativas/ iteracoes
SNRdB = [0:25]; %SNR em dB
SNR = 10.^(SNRdB/10);

c_SISO = zeros(1,length(SNR)); %capacidade do sistema SISO
c_SIMO = zeros(1,length(SNR)); %capacidade do sistema SIMO
c_MISO = zeros(1,length(SNR)); %capacidade do sistema MISO
c_MIMO = zeros(1,length(SNR)); %capacidade do sistema MIMO
c_MIMO4 = zeros(1,length(SNR));%capacidade do sistema MIMMO 4x4

for ite = 1:ITER
h_SISO = (randn +j*randn)/sqrt(2);              %
h_SIMO = (randn(N,1)+j*randn(N,1))/sqrt(2);     %
h_MISO = (randn(1,M)+j*randn(1,M))/sqrt(2);     %
h_MIMO = (randn(N,M)+j*randn(N,M))/sqrt(2);     %
h_MIMO = (randn(N,M)+j*randn(N,M))/sqrt(2);     %
for K = 1:length(SNR)
c_SISO(K) = c_SISO(K) + log2(1+ SNR(K)*norm(h_SISO)^2);
c_SIMO(K) = c_SIMO(K) + log2(1+ SNR(K)*norm(h_SIMO)^2);
c_MISO(K) = c_MISO(K) + log2(1+ SNR(K)*norm(h_MISO)^2/M);
c_MIMO(K) = c_MIMO(K) + log2(det(eye(N)+SNR(K)*h_MIMO*h_MIMO'/M));
end
end

c_SISO = c_SISO/ITER;
c_SIMO = c_SIMO/ITER;
c_MISO = c_MISO/ITER;
c_MIMO = c_MIMO/ITER;

plot(SNRdB,c_SISO,'r - .',SNRdB,c_SIMO,'b - o',SNRdB,c_MISO,'m',SNRdB,c_MIMO,'k - *')

legend('SISO 1x1','SIMO 1x2','MISO 2x1','MIMO 2x2')
xlabel('SNR em dB')
ylabel('Capacidade (b/s/Hz)')
title('Capacidade Vs. SNR')
grid;