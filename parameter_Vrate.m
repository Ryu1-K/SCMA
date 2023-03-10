function[rate] = parameter_Vrate(F, spa)

%%  Vの0の割合  % a=00, b=01, c=10, d=11
if spa < 0.5                                                %sparsity５０以下
    rate.a = 0;
    rate.b = floor(F * spa + 0.5);             %含まれる0の個数
    rate.c = rate.b;
    rate.d = F - (rate.b * 2);                                %他全部1
elseif spa == 0.5
    rate.a = 0;
    rate.b = floor(F * spa);
    rate.c = ceil(F * spa);
    rate.d = 0;
elseif spa > 0.5                                %50%以上
    spa_Im = 1-spa;
    rate.b = floor(F * spa_Im + 0.5);          %含まれる1の個数
    rate.c = rate.b;
    rate.a = F - (rate.b * 2);                  %他全部０
    rate.d = 0;
end