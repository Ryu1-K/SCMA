function[A, V, num, V_warn] = V_gen(rate, F,  L)

A_tmp = (1/2)*[-1 -1; -1 1];
A = sparse( kron(eye(L), A_tmp));

N = 2;              %�ʂ肩���@����N�ɂ��Ȃ�

num.V00 = rate.a;       %���Ӗ��Ȃ��Ƃ��Ă܂�
num.V01 = rate.b;
num.V10 = rate.c;
num.V11 = rate.d;

v00 = [0 0];
v01 = [0 1];
v10 = [1 0];
v11 = [1 1];

V00 = repmat(v00, num.V00, 1);
V01 = repmat(v01, num.V01, 1);
V10 = repmat(v10, num.V10, 1);
V11 = repmat(v11, num.V11, 1);

V_or = vertcat(V00, V01, V10, V11);

%%  V�����J�n
%%  ��点�Ȃ��ꍇ
V_warn = 0;  %V_l����������`����ϐ�
loop_ct = 0;
while(1)
    loop_ct = loop_ct + 1;
    
    if loop_ct == 2
        %             fprintf('�P��ڂ̔�蔭��\n');
        V_warn = 1;
    end
    
    if loop_ct == 100                                           %�����E�o
        fprintf('V_l�̔��Ȃ����͍̂��܂���\n');
        V_warn = 99;
        break
    end
    
    ran = randperm(F);
    va = 1:1:F;
    vb = F+1:1:2*F;
    vc = ran + F;
    V = zeros(F, N);
    V(va) = V_or(ran);
    V(vb) = V_or(vc);                                                %�ŏ���1�ߐ���
    V_l = cell(1, L);
    V_l{1} = V;
    
    exit = 1;                                                   %�E�o�ł���O��
    for V_ct = 2:L                                               %�Ίp�s��Ɋg�U
        ran = randperm(F);
        vc = ran + F;
        V_l{V_ct} = zeros(F, N);
        V_l{V_ct}(va) = V_or(ran);
        V_l{V_ct}(vb) = V_or(vc);
        V = sparse(blkdiag(V, V_l{V_ct}));
    end
    
    %% ����Ă�̂��Ȃ����`�F�b�N
        chk = 2;
        for chk_ct = 1:L-1
            for chk_ct2 = chk:L
                if V_l{chk_ct} == V_l{chk_ct2}
                    exit = 0;                                               %���������蒼��
                end
            end
            chk = chk + 1;
        end
     %%
    if exit == 1
        %             fprintf('���Ȃ�V_l��������\n');
        %             fprintf('V_l���������[�v��:%d\n', loop_ct);
        break;                                                              %���[�v�E�o
    end
end
