%% Script :SCMA-GaBP
clear
clc

%% Basic Simulation Setting
SIM.nworker = 2;        % Number of cores
SIM.det     = 1;         % Detector( 0:pda 1:gabp )
SIM.wloop   = 5;
SIM.nloop   = 10^SIM.wloop;     % Number of transmission
SIM.errmax  = SIM.nloop/10; 
SIM.K       = 16;        % Number of repetition
SIM.EsN0    = 10:1:24;    % SNR
SIM.EnSW    = 1;         % 0: Eb/N0, 1: Ex/N0, 
SIM.ndata   = 2;


%% Option Standard Setting
SIM.weight  = 1;    % Weight ( 0:MF-GaBP 1:LLR-GaBP)
SIM.damp    = 2;    % Damped ( 0:Conv 1:FIR 2:IIR)
SIM.method  = 0;    % Method ( 0:Conv 1:Normalizing 2:Limiting)
SIM.a       = 0;    % Parameter : Normalizing, Limiting
SIM.dya     = 0;
SIM.d       = 0.5;  % Parameter : Dampping factor

%% Simulation Script
SC.L = 72;      % Num. of layer (Num. of users)
SC.F = 48;      % Num. of dimension of sc
SC.M = 2;       % Num. of dimension of md
SC.N = 1;       % Num. of Rx antennas

%%  スパーシティ設定
SC.active = 0.3; 
SC.spa = (2-SC.active)/2;
SC.rate= parameter_Vrate(SC.F, SC.spa);  

%% Execute simulation 
% Num. of demensions
CH.nt = SC.M*SC.L;
CH.nr = SC.N*SC.F;

tic;
RES = zeros(SIM.K, length(SIM.EsN0),4);
% for idx_worker = 1:SIM.nworker %% debugのときはこっち
parfor idx_worker = 1:SIM.nworker
    RES_ = main_task(SIM, CH, SC, idx_worker);
    RES = RES +  RES_;
end
toc;
    
%% Result
SIM.FER = RES(:,:,1)./RES(:,:,2);
SIM.BER = RES(:,:,3)./RES(:,:,4);

plot_ber;
SIM.SC = SC;

fn = ['DATA\scma_' int2str(SC.L) '_' int2str(SC.F) '_' int2str(SC.M) '_' int2str(SC.active*10) '_' int2str(SIM.EnSW) '_' int2str(SIM.wloop)];
save(fn,'SIM');