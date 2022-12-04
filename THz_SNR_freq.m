%% params
clear all; clc;
B = 50e9; % 50 GHz
PSD=1e-22/(B); No=2*B*PSD;
d1 = [5,30,50]; % m
d2 = 5;
% d2 = [5, 30, 50];
Fdb = 10; % dB
F = 10^(Fdb/10);
T = 290; % K
k = 1.381e-23;
Ptx = 1; % W
N=64; 

k_rice = 5;

numSamples = 2000;
f_Thz = linspace(100, 1000, numSamples) * 10^9; % Hz
fc = f_Thz;
HITRANparams = importdata('data_freq_abscoe.txt');
loss_SR = zeros(numSamples, length(d1));
snr_los = zeros(numSamples, length(d1));
lossSpreadDb = zeros(numSamples, length(d1));
lossAbsDb = zeros(numSamples, length(d1));

%%
for freqIndex = 1:numSamples
	for distIndex = 1:length(d1)
%% parameters    
    [g1,g2]=chan_gen(fc(freqIndex),N,k_rice);
    alpha_i=abs(g1); theta_i = angle(g1);
    beta_i=abs(g2); psi_i = angle(g2);
    lambda=3e8/fc(freqIndex);
    d_IRS=0.3*lambda;
    zeta_PL=2*sqrt(pi*2)*d_IRS^2/(lambda^2); 
    zeta_PL_dB = mag2db(zeta_PL);
    
%% LoS loss
    lossSpreadDb(freqIndex, distIndex) = getSpreadLoss(f_Thz(freqIndex), sqrt(d1(distIndex)^2+d2^2));
	lossAbsDb(freqIndex, distIndex) = getAbsLoss(f_Thz(freqIndex), sqrt(d1(distIndex)^2+d2^2), HITRANparams);
	loss_los(freqIndex, distIndex) = lossSpreadDb(freqIndex, distIndex) + lossAbsDb(freqIndex, distIndex);

%% IRS loss
    lossSpreadDb(freqIndex, distIndex) = getSpreadLoss(f_Thz(freqIndex), d1(distIndex));
	[lossAbsDb(freqIndex, distIndex),kfParam_SR] = getAbsLoss(f_Thz(freqIndex), d1(distIndex), HITRANparams);
	loss_SR(freqIndex, distIndex) = lossSpreadDb(freqIndex, distIndex) + lossAbsDb(freqIndex, distIndex);
		
    lossSpreadDb2 = getSpreadLoss(f_Thz(freqIndex), d2);
	[lossAbsDb2,kfParam_RD] = getAbsLoss(f_Thz(freqIndex), d2, HITRANparams);
	loss_RD = lossSpreadDb2 + lossAbsDb2;
    
    k_thz(freqIndex,distIndex) = kfParam_RD;
    c = physconst('LightSpeed');
    l_SR = c/(4*pi*fc(freqIndex)*d1(distIndex))*exp(-0.5*kfParam_SR*d1(distIndex))*exp(-1j*2*pi*fc(freqIndex)*d1(distIndex));
    l_RD = c/(4*pi*fc(freqIndex)*d2)*exp(-0.5*kfParam_RD*d2)*exp(-1j*2*pi*fc(freqIndex)*d2);
   
%% SNR calculation
L4 = laguerreL(1/2,-k_rice)^4;
L_T = (abs(sum(zeta_PL*l_SR*l_RD)));
psi = L4*L_T^2*(1+k_rice^2);
snr_avg_rice(freqIndex,distIndex)=(N^2*pi^2*psi+L_T^2*N*(16*(1+k_rice^2)-pi^2*L4))...
    *Ptx/(16*2e-18);

snr_avg(freqIndex,distIndex)=(abs(sum(zeta_PL*l_SR*l_RD)))^2*...
        (N^2*pi^2+N*(16-pi^2))*Ptx...
        /(16*2e-18);
    
    
    end
end

snr_avg_db = 10*log10(snr_avg);
snr_rice_db = 10*log10(snr_avg_rice);

%% fig2

figure('DefaultAxesFontSize',18);
for distIndex = 1:length(d1)
	plot(f_Thz/1e9, (snr_rice_db(:, distIndex)),'-o',...
        'MarkerSize',10,...
        'MarkerIndices',1:100:length(f_Thz/1e9),'linewidth',4) 
    hold on
end

for distIndex = 1:length(d1)
	plot(f_Thz/1e9, (snr_avg_db(:, distIndex)),'-|',...
        'MarkerSize',10,...
        'MarkerIndices',1:100:length(f_Thz/1e9), 'linewidth',2) 
    hold on
end

yticks([-100 -75 -50 -25 0 25 50])
% ylim([-100 50])
xlim([fc(1)/1e9 fc(end)/1e9])
legend("d1 = " + d1(1) + " m, k = " + k_rice,"d1 = " + d1(2) + " m, k = " + k_rice,"d1 = " + d1(3) + " m, k = " + k_rice,...
    "d1 = " + d1(1) + " m","d1 = " + d1(2) + " m, IRS","d1 = " + d1(3) + " m, IRS");
ylim([-100 100])
hold on;
xlabel("Frequency (GHz)"); ylabel("SNR (dB)");
title("N="+N)
grid on
grid minor


%% k
% figure('DefaultAxesFontSize',18);
% for distIndex = 1:length(d1)
% 	plot(f_Thz/1e9, (k_thz(:, distIndex))/(100/3.4),'k-',...
%         'MarkerSize',10,...
%         'MarkerIndices',1:100:length(f_Thz/1e9),'linewidth',2)
%     hold on
% end
% xlabel("Frequency (GHz)"); ylabel("Absorption Coefficient");
% xlim([100 1000])

%% fig 3
% for distIndex = 1:length(d1)
% 	semilogy(f_Thz/1e9, (ber_thz1(:, distIndex)),'->',...
%         'MarkerSize',10,...
%         'MarkerIndices',1:100:length(f_Thz/1e9),'linewidth',2) 
%     hold on
% end
% 
% 
% for distIndex = 1:length(d1)
% 	semilogy(f_Thz/1e9, (ber_thz2(:, distIndex)),'-o',...
%         'MarkerSize',10,...
%         'MarkerIndices',1:100:length(f_Thz/1e9), 'linewidth',2) 
%     hold on
% end
% legend("d1 = " + d1(1) + " m, BPSK","d1 = " + d1(2) + " m, BPSK","d1 = " + d1(3) + " m, BPSK",...
%     "d1 = " + d1(1) + " m, QPSK","d1 = " + d1(2) + " m, QPSK","d1 = " + d1(3) + " m, QPSK");
% xlabel("Frequency (GHz)"); ylabel("BER");
% title("N="+N)
% grid on
% grid minor
