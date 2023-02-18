function[rate] = parameter_Vrate(F, spa)

%%  V‚Ì0‚ÌŠ„‡  % a=00, b=01, c=10, d=11
if spa < 0.5                                                %sparsity‚T‚OˆÈ‰º
    rate.a = 0;
    rate.b = floor(F * spa + 0.5);             %ŠÜ‚Ü‚ê‚é0‚ÌŒÂ”
    rate.c = rate.b;
    rate.d = F - (rate.b * 2);                                %‘¼‘S•”1
elseif spa == 0.5
    rate.a = 0;
    rate.b = floor(F * spa);
    rate.c = ceil(F * spa);
    rate.d = 0;
elseif spa > 0.5                                %50%ˆÈã
    spa_Im = 1-spa;
    rate.b = floor(F * spa_Im + 0.5);          %ŠÜ‚Ü‚ê‚é1‚ÌŒÂ”
    rate.c = rate.b;
    rate.a = F - (rate.b * 2);                  %‘¼‘S•”‚O
    rate.d = 0;
end