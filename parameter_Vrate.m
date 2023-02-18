function[rate] = parameter_Vrate(F, spa)

%%  V��0�̊���  % a=00, b=01, c=10, d=11
if spa < 0.5                                                %sparsity�T�O�ȉ�
    rate.a = 0;
    rate.b = floor(F * spa + 0.5);             %�܂܂��0�̌�
    rate.c = rate.b;
    rate.d = F - (rate.b * 2);                                %���S��1
elseif spa == 0.5
    rate.a = 0;
    rate.b = floor(F * spa);
    rate.c = ceil(F * spa);
    rate.d = 0;
elseif spa > 0.5                                %50%�ȏ�
    spa_Im = 1-spa;
    rate.b = floor(F * spa_Im + 0.5);          %�܂܂��1�̌�
    rate.c = rate.b;
    rate.a = F - (rate.b * 2);                  %���S���O
    rate.d = 0;
end