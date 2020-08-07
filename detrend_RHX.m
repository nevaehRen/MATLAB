

function [dt_ecgnl] = detrend_RHX(Ca_red,ShowOrnot)
% this function is used to remove the basal noise
% dt_ecgl linear denoise
% dt_ecgnl non-linear denoise
% and normlize trace to 0-1

ecgnl    = Ca_red;
t        = [1:length(ecgnl)]';

%% method 1 : detrend function
dt_ecgl  = detrend(ecgnl);

%% method 2 : polyfit remove basal function
% opol     = 6;
% [p,s,mu] = polyfit(t,ecgnl,opol);
% basal    = polyval(p,t,[],mu);

%% method 3 : exponenital decay remove basal function
f = fittype('a*exp(-x/b)+c');    
[fun_basal,    ~] = fit(t,ecgnl,   f,'StartPoint',[max(ecgnl) t(end)/2 min(ecgnl)],'Lower',[0 0 0],'Upper',[max(ecgnl)  t(end) max(ecgnl)]);
            basal = fun_basal.a.*exp(-[t]/fun_basal.b)+fun_basal.c;
% figure; hold on;
% plot(t, basal ,'g','linewidth',2)
% plot(t,Ca_red,'k','linewidth',2)

%%
dt_ecgnl = ecgnl - basal;
dt_ecgnl = (dt_ecgnl-min(dt_ecgnl))/(max(dt_ecgnl-min(dt_ecgnl)));

if nargin>1
    
figure;hold on;

subplot(3,1,1)
plot(t,dt_ecgl), grid;
title 'Detrended  Signals', ylabel 'Ca';

subplot(3,1,2)
plot(t,dt_ecgnl), grid;
xlabel Sample, ylabel 'Ca';


subplot(3,1,3);hold on;
plot(t,basal,'r'), grid;
plot(t,ecgnl,'k'), grid;
xlabel Sample, ylabel 'Ca';

end

