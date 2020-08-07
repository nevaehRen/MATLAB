
function Peaks = Peaks_with_Period_RHX(Data,Params,Name)


%% ===>> Step1: peaks searching parameters
if nargin<2
    Prominence = 0.2;
    Distance   = 2;
    MinPeakHeight = .2;
    Params.Prominence =Prominence;
    Params.Distance =Distance;
    Params.MinPeakHeight =MinPeakHeight;
else
    Prominence = Params.Prominence;
    Distance = Params.Distance;
    MinPeakHeight = Params.MinPeakHeight;
end

%% ===>> Step2: show peaks searching result, first row using halfheight, second row using halfprominence, third row using

figure;set(gcf,'Position',[50  200  1500 600], 'color',[1 1 1]);
ax(1) = subplot(3,1,1);
findpeaks(Data,'MinPeakProminence',Prominence,'MinPeakDistance',Distance,'MinPeakHeight',MinPeakHeight,'WidthReference','halfheight','Annotate','extents')
ax(2) = subplot(3,1,2); hold on; plot(Data,'k');
[pks,locs,w,p,bxPk,byPk,wxPk]      = findpeaks_RHX(Data,'MinPeakProminence',Prominence,'MinPeakDistance',Distance,'MinPeakHeight',MinPeakHeight,'WidthReference','halfheight','Annotate','extents');
%{  
 ====  peaks parameters ====
pks:  local maxima peak (y)
locs: peak location (x)
w:    width of the peaks
p:    prominences 
bxPk: local minimal (x, start and end)
byPk: local minimal (y)
wxPk: half width (x, start and end)
%}


%% ===>> Step3: recording peak features

Peaks(1).Params = Params;
Peaks(1).Data   = Data;

% --- Step3.1: adding peak(t3) and start(t1) and end(t5) of the local minimal 
for i=1:length(locs)
    % t1: minimal start point
    Peaks(i).t1 = bxPk(i,1); 
    Peaks(i).s1 = Data(Peaks(i).t1); 
    
    % t3: peak point
    Peaks(i).t3 = locs(i);
    Peaks(i).s3 = Data(Peaks(i).t3);
    
    % t5: minimal end point
    Peaks(i).t5 = bxPk(i,2);
    Peaks(i).s5 = Data(Peaks(i).t5);
end

% --- Step3.2: adding half increase(t2) and half decrease(t4), as well as 25% 50%
% 75% increase and decrease
Peaks = Fresh_Peaks_t2_t4_t75_t50_t25(Peaks);

% --- Step3.3: adding half increase(t2) and half decrease(t4), as well as 25% 50%
for i=1:length(locs)

    Peaks(i).t2_halfheight = wxPk(i,1);
    Peaks(i).s2_halfheight = interp1( Peaks(i).t1:Peaks(i).t3,Data(Peaks(i).t1:Peaks(i).t3) ,Peaks(i).t2_halfheight );
    
    Peaks(i).t4_halfheight = wxPk(i,2);
    Peaks(i).s4_halfheight = interp1( Peaks(i).t3:Peaks(i).t5,Data(Peaks(i).t3:Peaks(i).t5) ,Peaks(i).t4_halfheight );
    
    % adding all Time Difference 
    Peaks(i).t21 = Peaks(i).t2 - Peaks(i).t1;
    Peaks(i).s21 = Peaks(i).s2 - Peaks(i).s1;
    
    Peaks(i).t31 = Peaks(i).t3 - Peaks(i).t1;
    Peaks(i).s31 = Peaks(i).s3 - Peaks(i).s1;
    
    Peaks(i).t41 = Peaks(i).t4 - Peaks(i).t1;
    Peaks(i).s41 = Peaks(i).s4 - Peaks(i).s1;
    
    Peaks(i).t51 = Peaks(i).t5 - Peaks(i).t1;
    Peaks(i).s51 = Peaks(i).s5 - Peaks(i).s1;
    
    Peaks(i).t32 = Peaks(i).t3 - Peaks(i).t2;
    Peaks(i).s32 = Peaks(i).s3 - Peaks(i).s2;
    
    Peaks(i).t42 = Peaks(i).t4 - Peaks(i).t2;
    Peaks(i).s42 = Peaks(i).s4 - Peaks(i).s2;
    
    Peaks(i).t52 = Peaks(i).t5 - Peaks(i).t2;
    Peaks(i).s52 = Peaks(i).s5 - Peaks(i).s2;
    
    Peaks(i).t43 = Peaks(i).t4 - Peaks(i).t3;
    Peaks(i).s43 = Peaks(i).s4 - Peaks(i).s3;
    
    Peaks(i).t53 = Peaks(i).t5 - Peaks(i).t3;
    Peaks(i).s53 = Peaks(i).s5 - Peaks(i).s3;
    
    Peaks(i).t54 = Peaks(i).t5 - Peaks(i).t4;
    Peaks(i).s54 = Peaks(i).s5 - Peaks(i).s4;
    
    % adding Area 
    Peaks(i).A0 = sum(Data(Peaks(i).t1:Peaks(i).t5));
    Peaks(i).A1 = sum(Data(Peaks(i).t1:Peaks(i).t2));
    Peaks(i).A2 = sum(Data(Peaks(i).t2:Peaks(i).t3));
    Peaks(i).A3 = sum(Data(Peaks(i).t3:Peaks(i).t4));
    Peaks(i).A4 = sum(Data(Peaks(i).t4:Peaks(i).t5));  
    
    % adding increase time constant and delay time constant 
    Peaks(i).tau_up    = fit_Up_Peaks(Peaks,i);
    Peaks(i).tau_decay = fit_decay_Peaks(Peaks,i);
    % Peaks(i).tau_up = fit_Up_Peaks(Peaks,i,'show');
    % Peaks(i).tau_decay = fit_decay_Peaks(Peaks,i,'show');
    
    % adding peak shape features
    Peaks(i).oscillation_shape_area = sum(Data(Peaks(i).t1:Peaks(i).t3))/sum(Data(Peaks(i).t1:Peaks(i).t5));
    Peaks(i).oscillation_shape_T    = Peaks(i).tau_up/Peaks(i).tau_decay;
 
    % show and check 
    plot(Peaks(i).t1,Peaks(i).s1,'r^','markerfacecolor','r')
    plot(Peaks(i).t2,Peaks(i).s2,'r>','markerfacecolor','r')
    plot(Peaks(i).t3,Peaks(i).s3,'ro','markerfacecolor','r')
    plot(Peaks(i).t4,Peaks(i).s4,'r<','markerfacecolor','r')
    plot(Peaks(i).t5,Peaks(i).s5,'r^','markerfacecolor','r')
    
end


% --- Step3.4: adding basal Ca activity
if length(locs)~=0   
    T1 = [Peaks.t1];
    S1 = [Peaks.s1];
    T5 = [Peaks.t5];
    S5 = [Peaks.s5];
    Peaks(1).basal              = interp1([T1 T5(end)],[S1 S5(end)],1:length(Data))';
    Peaks(1).basal(1:T1(1))     = S1(1);
    Peaks(1).basal(T5(end):end) = S5(end);
else
    Peaks(1).basal = Data;
    Peaks.empty    = [];
end  


% --- Step3.5: adding period(defined by t2 time difference !)
if length(locs)>2
    for i=1:length(locs)-1
        Peaks(i).Period = Peaks(i+1).t2-Peaks(i).t2;
        Peaks(i).Fat    = Peaks(i).t42/Peaks(i).Period;
    end
end


% --- Step3.6: check, if has Name input, then save pics
ax(3) = subplot(3,1,3); hold on; plot(1:length(Data), Data,'k');
scatter([Peaks.T_up_25],[Peaks.s_up_25],[40], 'r>', 'filled');
scatter([Peaks.T_up_50],[Peaks.s_up_50],[40], 'g>', 'filled');
scatter([Peaks.T_up_75],[Peaks.s_up_75],[40], 'b>', 'filled');
scatter([Peaks.T_down_25],[Peaks.s_down_25],[40], 'r<', 'filled');
scatter([Peaks.T_down_50],[Peaks.s_down_50],[40], 'g<', 'filled');
scatter([Peaks.T_down_75],[Peaks.s_down_75],[40], 'b<', 'filled');
scatter([Peaks.t3],[Peaks.s3],[100], 'ko', 'filled');
scatter([Peaks.t2_halfheight],[Peaks.s2_halfheight],[40], 'mo');
scatter([Peaks.t4_halfheight],[Peaks.s4_halfheight],[40], 'mo');

linkaxes (ax,'x') ;

if nargin==3
    mkdir('Peaks');
    saveas(gcf,['Peaks/' Name '.png'])
    close();
else
    pause(10);
    close();
end



end
% clearvars -except Peaks Data Islet
% save([Islet(1).Name	'_Peaks.mat']);


function Peaks = Fresh_Peaks_t2_t4_t75_t50_t25(Peaks)
% fill t75 t50 and t25 for Peaks
% 25
% 50
% 75 point, all are based on the local maxmial and minmal

Data = Peaks(1).Data;

%% Step 1: Prepare the decreasing Ca data. interpolation help us find an accurate 25% 50% 75 decreasing points

for i=1:length(Peaks)   
    % Step 1.1 == using decresing Ca data
    T_ori  = (Peaks(i).t3:Peaks(i).t5);
    Ca     = Data(T_ori);
    T      = T_ori-T_ori(1);
    % Step 1.2 == 10x interpolation
    Amplitude_decrease = Ca(1)-Ca(end);
    T_detail           = 0:0.01:T(end);
    Ca_10x             = interp1(T,Ca,T_detail); 
    T_beta_75 = min(T_detail(Ca_10x <= (Amplitude_decrease*0.75+Ca(end)) ));
    T_beta_50 = min(T_detail(Ca_10x <= (Amplitude_decrease*0.5 +Ca(end)) ));
    T_beta_25 = min(T_detail(Ca_10x <= (Amplitude_decrease*0.25+Ca(end)) ));
    % Step 1.3 == real time 
    T_beta_75 = T_beta_75+T_ori(1);
    T_beta_50 = T_beta_50+T_ori(1);
    T_beta_25 = T_beta_25+T_ori(1);
    % Step 1.4 == record all decreasing 25% 50% 75 points
    Peaks(i).T_down_75= T_beta_75;   
    Peaks(i).T_down_50= T_beta_50;  
    Peaks(i).T_down_25= T_beta_25;  
    Peaks(i).s_down_75= (Amplitude_decrease*0.75+Ca(end));   
    Peaks(i).s_down_50= (Amplitude_decrease*0.5+Ca(end));  
    Peaks(i).s_down_25= (Amplitude_decrease*0.25+Ca(end));  
%{
     %%  === plot and check ===         
     figure;set(gcf,'Position',[50  200  1000 350], 'color',[1 1 1]);hold on;
     plot(T_ori,Ca,'b','linewidth',3)
     plot(T_beta_75,(Amplitude_decrease*0.75+Ca(end)),'r.','markersize',15,'linewidth',2)
     plot(T_beta_50,(Amplitude_decrease*0.5 +Ca(end)),'k.','markersize',15,'linewidth',2)
     plot(T_beta_25,(Amplitude_decrease*0.25+Ca(end)),'g.','markersize',15,'linewidth',2)
     set(gca,'linewidth',1.5 , 'Fontsize', 20, 'Fontname' , 'Times New Roman');
     plot(Peaks(i).t3,Peaks(i).s3,'bx','markersize',15,'linewidth',2)
     plot(Peaks(i).t4,Peaks(i).s4,'rx','markersize',15,'linewidth',2)
%}
end

%% Step 2: Prepare the increasing Ca data. interpolation help us find an accurate 25% 50% 75 increasing points

for i=1:length(Peaks)
    % Step 2.1 == using increasing Ca data
    T_ori  = (Peaks(i).t1:Peaks(i).t3);
    Ca     = Data(T_ori);
    T      = T_ori-T_ori(1);
    % Step 2.2 == 10x interpolation
    Amplitude_decrease = Ca(end)-Ca(1);
    T_detail = 0:0.01:T(end);
    Ca_10x = interp1(T,Ca,T_detail);
    T_beta_75=min(T_detail(Ca_10x >= (Amplitude_decrease*0.75+Ca(1)) ));
    T_beta_50=min(T_detail(Ca_10x >= (Amplitude_decrease*0.5+Ca(1))  ));
    T_beta_25=min(T_detail(Ca_10x >= (Amplitude_decrease*0.25+Ca(1)) ));
    % Step 2.3 == real time 
    T_beta_75 = T_beta_75+T_ori(1);
    T_beta_50 = T_beta_50+T_ori(1);
    T_beta_25 = T_beta_25+T_ori(1);
    % Step 2.4 == record all increasing 25% 50% 75 points
    Peaks(i).T_up_75= T_beta_75;   
    Peaks(i).T_up_50= T_beta_50;  
    Peaks(i).T_up_25= T_beta_25;  
    Peaks(i).s_up_75= (Amplitude_decrease*0.75+Ca(1));   
    Peaks(i).s_up_50= (Amplitude_decrease*0.5+Ca(1));  
    Peaks(i).s_up_25= (Amplitude_decrease*0.25+Ca(1)); 
%{
    % step3. according to expermintal f(x), calculate revised expermintal f_exper(x)         
    %          figure;set(gcf,'Position',[50  200  1000 350], 'color',[1 1 1]);hold on;
    %          plot(T_ori,Ca,'b','linewidth',3)
    %          plot(T_beta_75,(Amplitude_decrease*0.75+Ca(1)),'r.','markersize',15,'linewidth',2)
    %          plot(T_beta_50,(Amplitude_decrease*0.5 +Ca(1)),'k.','markersize',15,'linewidth',2)
    %          plot(T_beta_25,(Amplitude_decrease*0.25+Ca(1)),'g.','markersize',15,'linewidth',2)
    %          set(gca,'linewidth',1.5 , 'Fontsize', 20, 'Fontname' , 'Times New Roman');
    %          plot(Peaks(i).t3,Peaks(i).s3,'bx','markersize',15,'linewidth',2)
    %          plot(Peaks(i).t2,Peaks(i).s2,'rx','markersize',15,'linewidth',2)
%}
end

%% Step 3: t2 and t4 as half local minimal and maximal

for i=1:length(Peaks)
    Peaks(i).t2   = Peaks(i).T_up_50;  
    Peaks(i).s2   = Peaks(i).s_up_50; 
    Peaks(i).t4   = Peaks(i).T_down_50;  
    Peaks(i).s4   = Peaks(i).s_down_50; 
end
end




function tau = fit_decay_Peaks(Peaks,i,ShoworNot)
% this function fit decay time constant, using equation 'a*exp(-x/b)+c',

if i<=length(Peaks)
    Ca      = Peaks(1).Data;  % get Ca signal
    mini_T  = [Peaks(i).t3:Peaks(i).t5];  % using peak t3 to local minimal t5

    % fitting and get time constant tau
    if (Peaks(i).t5 - Peaks(i).t3)>1
        mini_T_shift = [mini_T-min(mini_T)];
        mini_Ca      = Ca(mini_T);
        f            = fittype('a*exp(-x/b)+c');    
        [cfun_WT,    ~] = fit(mini_T_shift(:),mini_Ca(:),   f,'StartPoint',[mini_Ca(1) mini_T_shift(end)/2 mini_Ca(end)],'Lower',[mini_Ca(end) 0 0],'Upper',[mini_Ca(1) mini_T_shift(end) mini_Ca(1)]);
        tau = cfun_WT.b;
    else
        tau = 0.5;
    end

    % check and show
    if nargin>2 & tau~=0.5
        S3      = Peaks(i).s3;
        S4      = interp1(mini_T,mini_Ca,Peaks(i).t4);
        figure;set(gcf,'Position',[50  200  1000 350], 'color',[1 1 1]);hold on;
        plot(Ca,'r','linewidth',2);
        plot(mini_T,interp1(mini_T,mini_Ca,mini_T),'b','linewidth',3)
        plot(Peaks(i).t3,S3,'b.','markersize',25)
        plot(Peaks(i).t4,S4,'r.','markersize',25)
        plot(mini_T,cfun_WT.a.*exp(-[mini_T-min(mini_T)]/cfun_WT.b)+cfun_WT.c,'g','linewidth',2)
        %xlim([50 150]);     
        set(gca,'linewidth',1.5 , 'Fontsize', 20, 'Fontname' , 'Times New Roman');
    end
end

end

function tau = fit_Up_Peaks(Peaks,i,ShoworNot)

% it gaves the increasing time constants
% using equation 'a*exp(-x/b)+c'

if i<=length(Peaks)
         Ca = Peaks(1).Data;
         mini_T  = [Peaks(i).t1:Peaks(i).t3]; % using peak t1 to local minimal t3
 
         if (Peaks(i).t3 - Peaks(i).t1)>1
             mini_T_shift = [mini_T-min(mini_T)];
             mini_Ca = Ca(mini_T(end:-1:1));
             f = fittype('a*exp(-x/b)+c');    
             [cfun_WT,    ~] = fit(mini_T_shift(:),mini_Ca(:),   f,'StartPoint',[mini_Ca(1) mini_T_shift(end)/2 mini_Ca(end)],'Lower',[mini_Ca(end) 0 0],'Upper',[mini_Ca(1) mini_T_shift(end) mini_Ca(1)]);
             tau = cfun_WT.b;
         else
             tau = 0.5;
         end

         if nargin>2 & tau~=0.5
             S3      = Peaks(i).s3;
             S4      = interp1(mini_T,mini_Ca,Peaks(i).t4);

             figure;set(gcf,'Position',[50  200  1000 350], 'color',[1 1 1]);hold on;
             plot(Ca,'r','linewidth',2);
             plot(mini_T(end:-1:1),interp1(mini_T,mini_Ca,mini_T),'b','linewidth',3)
             plot(Peaks(i).t3,S3,'b.','markersize',25)
             plot(Peaks(i).t4,S4,'r.','markersize',25)
             plot(mini_T,cfun_WT.a.*exp(-[mini_T(end:-1:1)-min(mini_T)]/cfun_WT.b)+cfun_WT.c,'g','linewidth',2)
             % xlim([50 150]);     
             set(gca,'linewidth',1.5 , 'Fontsize', 20, 'Fontname' , 'Times New Roman');
         end
end


end



