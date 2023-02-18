%% Plot BER (E_s/N0)

% SIM   : simulation parameter
%       SIM.nsamp : Number of transmissions
%       SIM.modu  : Select modulation 
%       SIM.K     : Number of repetition
%       SIM.dcL   : Select the way of ƒ¿ combination
%       SIM.EsN0  : SNR array

%% Script

FN = 'Times New Roman';
FS = 18;

h = semilogy(SIM.EsN0, SIM.BER(SIM.K,:));
% h = semilogy(SIM.EsN0, SIM.BER(2,:));
hold on

% axis([min(SIM.EsN0) max(SIM.EsN0) 10^(-6) 1]);
axis([min(SIM.EsN0) max(SIM.EsN0) 10^(-5) 1]);
grid on;
 
lb(1) = xlabel('\it{E}_{\rm{s}}\rm{/}\it{N}_{\rm{0}} \rm{[dB]}');
lb(2) = ylabel('BER');

% set(gca,'Linewidth',3.0,'FontName',FN,'FontSize',FS,'xTick',[min(SIM.EsN0):SIM.EsN0(2)-SIM.EsN0(1):max(SIM.EsN0)],'yTick',10.^[-6:0]);
set(gca,'Linewidth',3.0,'FontName',FN,'FontSize',FS,'xTick',[min(SIM.EsN0):4:max(SIM.EsN0)],'yTick',10.^[-6:0]);
set(h,'Linewidth',3.0);
set(lb,'Linewidth',3.0,'FontName',FN,'FontSize',FS);