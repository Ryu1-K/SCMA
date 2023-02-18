%% Post combine turbo equalization(QPSK)

% CH    : MIMO channel model
%       CH.nt : Number of transmit antennas
%       CH.nr : Number of receive antennas

% RX    : Receiver infomation

% SIM   : Simulation parameter
%       SIM.nsamp : Number of transmissions
%       SIM.modu  : Select modulation
%       SIM.K     : Number of repetition
%       SIM.dpL   : Select the way of create xh
%       SIM.dcL   : Select the way of ƒ¿ combination
%       SIM.EsN0  : SNR array

% SIM   : Option parameter
%       SIM.r1    : first refresh
%       SIM.r2    : refresh interval

%% Post combine turbo equalization

function DTC = det_gabp(CH, RX, SIM)

% Posterior LLR container
DTC.cllr = zeros(CH.nt, size(RX.y,2), SIM.K);

for idx_block = 1:size(RX.y,2)
    
    %%%%%%%%%% Basic Preparation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %     % Unweighted Average
    %     HH_abs = ( CH.H.*conj(CH.H) ).';
    
    %%%%%%%%%%% Detecton %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Symbol Replica matrix
    SR_matrix = zeros(CH.nt, CH.nr);
    
    for k = 1:SIM.K
        
        %%%%%%%%%% Calculate ƒ¿ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Cancell
        Canceller = sum( (CH.H.').*SR_matrix );
        yn = (RX.y(:,idx_block)).' - Canceller;
        Cancelled_matrix = yn(ones(CH.nt, 1),:);
        
        % Reconstruct matrix
        Reconstruct_matrix = (CH.H.').*SR_matrix;
        % Reconstruct and Mached Filtering
        hynm_matrix = (CH.H').*(Cancelled_matrix + Reconstruct_matrix);
        
        
        % Create First Weight
        if SIM.weight == 0               % MF-GaBP
            W = 2*sqrt(CH.Es)./CH.N0;
            
        elseif SIM.weight == 1           % LLR-GaBP
            dnm = CH.Es - SR_matrix.*conj(SR_matrix);
            element = CH.HH_abs.*dnm;
            hdh= sum(element);
            X = hdh(ones(CH.nt, 1), :) - element;
            W = 2*sqrt(CH.Es)./(CH.N0 + X);
        end
        
        % ƒ¿ matrix
        Alpha_matrix_c = W.*hynm_matrix;
        Alpha_ave_c = W.*CH.HH_abs;
        
        % Create final detection data
        DTC.cllr(:,idx_block,k) = sum(Alpha_matrix_c, 2);
        
        % ƒÀ matrix
        Beta_matrix_c = repmat( sum(Alpha_matrix_c, 2), [1, CH.nr]) - Alpha_matrix_c;
        Beta_ave_c = repmat( sum(Alpha_ave_c, 2), [1,CH.nr] )-Alpha_ave_c;
        
        % ƒÉ matrix
        % (1) Conv or Damp
        
        if SIM.damp == 0         % Conventional
            Damped_matrix = Beta_matrix_c;
            norm_num = Beta_ave_c;
            
        elseif SIM.damp == 1     % FIR Damping with damping factor
            if k < 3
                Damped_matrix = Beta_matrix_c;
                norm_num = Beta_ave_c;
            else
                Damped_matrix = SIM.d*Beta_matrix_c + (1-SIM.d)*Beta_matrix_p;
                norm_num = SIM.d*Beta_ave_c + (1-SIM.d)*Beta_ave_p;
            end
            
        elseif SIM.damp == 2     % IIR Damping with damping factor
            if k < 3
                Damped_matrix = Beta_matrix_c;
                norm_num = Beta_ave_c;
            else
                Damped_matrix = SIM.d*Beta_matrix_c + (1-SIM.d)*Damped_matrix;
                norm_num = SIM.d*Beta_ave_c + (1-SIM.d)*norm_num;
            end
            
        else
            disp('error: no implementation');
            exit
        end
        
        if SIM.dya == 1
            SIM.a = 1 + (1/SIM.K)*k;
        end
        
        
        % (2) Conv of Normalization or Limitation
        if SIM.method == 0
            Lambda_matrix = Damped_matrix;
            
        elseif SIM.method == 1  % Normalization
            % check 0
            norm_num(norm_num==0) = 1;
            % Normalization
            comp_ave = SIM.a./(sqrt(CH.Es)*norm_num);
            Lambda_matrix = comp_ave.*Damped_matrix;
            
        elseif SIM.method == 2  % Limitation
            Damped_matrix_r = real(Damped_matrix);
            Damped_matrix_i = imag(Damped_matrix);
            
            Damped_matrix_r(Damped_matrix_r > SIM.a) =  SIM.a;
            Damped_matrix_r(Damped_matrix_r <-SIM.a) = -SIM.a;
            
            Damped_matrix_i(Damped_matrix_i > SIM.a) =  SIM.a;
            Damped_matrix_i(Damped_matrix_i <-SIM.a) = -SIM.a;
            
            Lambda_matrix = Damped_matrix_r + 1j*Damped_matrix_i;
            
        elseif SIM.method == 3  % staticaly scaling
            Lambda_matrix = (1/SIM.a)*Damped_matrix;
            
        else
            disp('error: no implementation');
            exit
        end
        
        % Create Symbol Replica
        SR_matrix = tanh( real(Lambda_matrix))*sqrt(CH.Es);
        
        % Update
        Beta_matrix_p = Beta_matrix_c;
        Beta_ave_p = Beta_ave_c;
        
    end
end

end
