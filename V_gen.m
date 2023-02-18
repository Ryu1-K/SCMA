function[A, V, num, V_warn] = V_gen(rate, F,  L)

A_tmp = (1/2)*[-1 -1; -1 1];
A = sparse( kron(eye(L), A_tmp));

N = 2;              %ぬりかえ　実はNによらない

num.V00 = rate.a;       %無意味なことしてます
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

%%  V生成開始
%%  被らせない場合
V_warn = 0;  %V_l生成成功を伝える変数
loop_ct = 0;
while(1)
    loop_ct = loop_ct + 1;
    
    if loop_ct == 2
        %             fprintf('１回目の被り発生\n');
        V_warn = 1;
    end
    
    if loop_ct == 100                                           %強制脱出
        fprintf('V_lの被らないものは作れません\n');
        V_warn = 99;
        break
    end
    
    ran = randperm(F);
    va = 1:1:F;
    vb = F+1:1:2*F;
    vc = ran + F;
    V = zeros(F, N);
    V(va) = V_or(ran);
    V(vb) = V_or(vc);                                                %最初の1個め生成
    V_l = cell(1, L);
    V_l{1} = V;
    
    exit = 1;                                                   %脱出できる前提
    for V_ct = 2:L                                               %対角行列に拡散
        ran = randperm(F);
        vc = ran + F;
        V_l{V_ct} = zeros(F, N);
        V_l{V_ct}(va) = V_or(ran);
        V_l{V_ct}(vb) = V_or(vc);
        V = sparse(blkdiag(V, V_l{V_ct}));
    end
    
    %% 被ってるのがないかチェック
        chk = 2;
        for chk_ct = 1:L-1
            for chk_ct2 = chk:L
                if V_l{chk_ct} == V_l{chk_ct2}
                    exit = 0;                                               %あったらやり直し
                end
            end
            chk = chk + 1;
        end
     %%
    if exit == 1
        %             fprintf('被らないV_l生成成功\n');
        %             fprintf('V_l生成時ループ回数:%d\n', loop_ct);
        break;                                                              %ループ脱出
    end
end
