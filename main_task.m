function RES = main_task(SIM, CH, SC, idx_worker)
% Create random stream( by Mersenne twister ) and set global stream
s = RandStream('mt19937ar','Seed', idx_worker);
RandStream.setGlobalStream(s);
% Energy per Symbol (equivalent to bit in SCMA)
CH.Es = 1;
% Variance of Fading
CH.Nh = 1;
% Index of BER array
sim_ct = 1;

for En = SIM.EsN0
    ERR.noe = zeros(SC.L,SIM.K);        ERR.nod = zeros(SC.L,SIM.K);
    ERR.noef = zeros(SC.L,SIM.K);       ERR.nodf = zeros(SC.L,SIM.K);
    
    % Calculate N0 from SNR (Variance of Noise)
    switch(SIM.EnSW)
        case 0
            CH.N0 = 10^(-En/10)*CH.Es*((SC.rate.b+SC.rate.c+2*SC.rate.d)*SC.L/2)/(SC.L*SC.M);
        case 1
            CH.N0 = 10^(-En/10)*CH.Es*((SC.rate.b+SC.rate.c+2*SC.rate.d)*SC.L/2)/SC.L;
    end
    for idx_loop = 1:ceil(SIM.nloop/SIM.nworker)
        [SC.A, SC.V, num, V_warn] = V_gen(SC.rate, SC.F, SC.L);
        % Create H matrix
        CH.H = create_H_matrix(CH,SC);
        % Unweighted Average
        CH.HH_abs = ( CH.H.*conj(CH.H) ).';
        
        %%%%%%%%%% Transmitter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        TX.b = randi([0 1], SC.L,SIM.ndata);
        TX.bb = zeros(SC.L*SC.M, SIM.ndata/SC.M);
        TX.bb(1:2:end,:)=TX.b(:,1:2:end);
        TX.bb(2:2:end,:)=TX.b(:,2:2:end);
        TX.x = 2*TX.bb-1;
        
        %%%%%%%%%% Channel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        RX.z = ( randn(CH.nr, size(TX.x,2)) + 1i*randn(CH.nr,size(TX.x,2)) )*sqrt(CH.N0/2);
        RX.y = CH.H*TX.x + RX.z;
        
        %%%%%%%%%% Receiver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        DTC = det_gabp(CH, RX, SIM);
        
        RX.cllr = zeros(SC.L, SIM.ndata,SIM.K);
        RX.cllr(:,1:2:end,:) = DTC.cllr(1:2:end,:,:);
        RX.cllr(:,2:2:end,:) = DTC.cllr(2:2:end,:,:);
        
        % Hard Dec.
        RX.b  = ( RX.cllr>0 );
        tmp = reshape(sum(RX.b ~= TX.b,2),SC.L,SIM.K);
        ERR.noe = ERR.noe + reshape(sum(RX.b ~= TX.b,2),SC.L,SIM.K);
        ERR.nod = ERR.nod + SIM.ndata;
        ERR.noef = ERR.noef + (tmp~=0);
        ERR.nodf = ERR.nodf + 1;
        if(mean(ERR.noef(:,end))>(SIM.errmax/SIM.nworker))
            break
        end
    end
    % Store FER
    RES(:,sim_ct,1) = sum(ERR.noef)';       RES(:,sim_ct,2) = sum(ERR.nodf)';
    % Store BER (each llr)
    RES(:,sim_ct,3) = sum(ERR.noe)';        RES(:,sim_ct,4) = sum(ERR.nod)';
    sim_ct = sim_ct + 1;
end
end