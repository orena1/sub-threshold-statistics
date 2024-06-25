
clear
% Load data
addpath('/m_files')
yaml_file='par.yaml';
y_str = ReadYaml(yaml_file);

   
ne=y_str.ne;ni=y_str.ni;
dt=y_str.dt;
TstopStim= y_str.time_sim;
Tstart=y_str.Tstart; 
[m1St]=Tstart/dt; [m1En]=(TstopStim-1)/dt;
time_sim =y_str.time_sim;

cd Data
% cd DataFigurePaper
iter=1;
[x,y]=textread(['spikee' num2str(iter) '.txt'],'%f %f');[x2,y2]=textread(['spikei' num2str(iter) '.txt'],'%f %f');
% [x,y]=textread(['spikee.txt'],'%f %f');[x2,y2]=textread(['spikei.txt'],'%f %f');

figure;plot(x(1:length(y)),y(1:length(y)),'.r' );hold on;plot(x2(1:length(y2)),y2(1:length(y2)),'.' );
% Rates Ver1
Se=[x(1:length(y)),y(1:length(y))];
Si=[x2(1:length(y2)),y2(1:length(y2))];
Si(:,2)=Si(:,2)-ne+1;
% Se(:,1)=Se(:,1)+1;
Se(:,2)=Se(:,2)+1;

Se_in=[];
Se_in(:,1)=round(Se(:,1)./dt)+1;
Se_in(:,2)=round(Se(:,2));
Se_out=sparse(Se_in(:,2),Se_in(:,1),ones(length(Se_in),1));

Si_in(:,1)=round(Si(:,1)./dt)+1;
Si_in(:,2)=round(Si(:,2));
Si_out=sparse(Si_in(:,2),Si_in(:,1),ones(length(Si_in),1));


%Stim 1
Si_1=Si_out(:,m1St:m1En);
FrI=sum(Si_1,2)./(TstopStim-Tstart)*100;% num of spikes in ms.. 1 tau
mean(FrI)
Se_1=Se_out(:,m1St:m1En);
FrE=sum(Se_1,2)./(TstopStim-Tstart)*100;% num of spikes in ms.. 1 tau
results.mFrE=mean(FrE);results.mFrI=mean(FrI);
results.sFrE=std(FrE);results.sFrI=std(FrI);


% disp(['FrE:' num2str(mean(FrE)),'+-',num2str(std(FrE)) '  FrI:' num2str(mean(FrI)),'+-',num2str(std(FrI))])

% Calc CV
[results.mCVE,results.sCVE, results.mCVI,results.sCVI,CVEvec,CVIvec] = CV_moments(Se_1,Si_1);
% disp(['CVE:' num2str(mCVE),'+-',num2str(sCVE) 'CVI:' num2str(mCVI),'+-',num2str(sCVI)])

%%

tte=(1:size(Se_out,2))*dt;
tti=(1:size(Si_out,2))*dt;
figure;plot(tte,mean(Se_out)./dt*100);hold all;plot(tti,mean(Si_out)./dt*100)

%% Voltage and rate Ver2

WBFlag=0;
CutSt=0;
ve=load('ve1.txt');
[results.mVE,results.sdVE, results.mSDVE,results.sdSDVE, results.SKVE, results.V1E, results.sdV1E, results.skV1E, FrEfromV, VEpop,sdEpop] = voltage_moments(ve,WBFlag,m1St,m1En,CutSt);
vi=load('vi1.txt');
vi(:,ni+1:end)=[];
[results.mVI,results.sdVI, results.mSDVI,results.sdSDVI, results.SKVI, results.V1I, results.sdV1I, results.skV1I, FrIfromV,VIpop,sdIpop] = voltage_moments(vi,WBFlag,m1St,m1En,CutSt);
mean(FrEfromV)*100
mean(FrIfromV)*100
ve(m1En:end,:)=[];
vi(m1En:end,:)=[];
if CutSt
ve(1:m1St,:)=[];
vi(1:m1St,:)=[];
end

 V=ve(:,2:end);V(1:10,:)=[]; V(logical(V>-40))=V(find(V>-40)-1);mean(skewness(V))
%% Fit single neuron to a Gaussian distribution
ii=8
figure
thresh=-40;
  T=ve(:,1);V=ve(:,ii);  V(1)=[];T(1)=[];
  V(logical(V>thresh))=V(find(V>thresh)-1);
  V(1:end/10)=[];
 sk=skewness(V); 

x=V;
pd=fitdist(x,'Normal');
BIN_MAX=-20;BIN_WIDTH=0.5;
BIN_RANGE=-80:BIN_WIDTH:BIN_MAX;
% Empirical distribution
y_data=hist(x,BIN_RANGE);
% bar(BIN_RANGE,y_data./sum(y_data)/BIN_WIDTH,'facealpha',.5,'edgecolor','none')

bar(BIN_RANGE,y_data./sum(y_data)/BIN_WIDTH)

hold on
% x_values = BIN_RANGE;
x_values = -90:0.01:BIN_MAX;
y = pdf(pd,x_values);
plot(x_values,y,'LineWidth',1)
hold on
plot([thresh thresh],[ 0 0.2],'-r')
xlim([-80 -30])
axis square; box off;
title(num2str(sk))


%%
figure
subplot(2,1,1)
% plot(ve(:,1),ve(:,2:end)-repmat(ve(39/dt,2:end),size(ve,1),1))

plot(ve(:,1),ve(:,2:end))
hold all
plot(ve(1:length(VEpop),1),VEpop,'k','LineWidth',4)
plot(ve(1:length(sdEpop),1),sdEpop,'g','LineWidth',4)
subplot(2,1,2)
plot(vi(:,1),vi(:,2:40))
hold all
plot(vi(:,1),VIpop,'k','LineWidth',4)
plot(vi(:,1),sdIpop,'g','LineWidth',4)
%% Plot log-normal and CV Vs Fr
% addpath(genpath('/Users/darshanr/data'))

figure(1)
plot(results.V1E,results.sdV1E,'o')
hold on
xlim([-70 -30]);ylim([0 6]);
axis square;box off;
x=results.V1E';
y=results.sdV1E';
[Rsq,yCalc]= MyLinearReg(x,y); 
Rsq
[r, pr]=corrcoef(x,y)
plot(x,yCalc,'r')
title(sprintf('mean voltage Exc: %.1f +- %.1f; SD= %.1f',results.mVE,results.sdVE,results.mSDVE))
ylim([0 5])

figure(2)
plot(results.V1I,results.sdV1I,'o')
hold on
xlim([-70 -30]);ylim([0 6]);
axis square;box off;
x=results.V1I';
y=results.sdV1I';
plot(x,y,'o')

[Rsq,yCalc]= MyLinearReg(x,y);
Rsq
[r, pr]=corrcoef(x,y)
plot(x,yCalc,'r')

title(sprintf('mean voltage Inh: %.1f +- %.1f; SD= %.1f',results.mVI,results.sdVI,results.mSDVI))
ylim([0 5])

% Empirical distribution
figure(5)
% BIN_WIDTH=.5;
BIN_WIDTH=1;
BIN_MAX=-30;
BIN_RANGE=[-70:BIN_WIDTH:-30];
x=results.V1E;


BIN_RANGE = -80:BIN_WIDTH:BIN_MAX;
% x_values = -90:0.01:BIN_MAX;
y_data=hist((x),BIN_RANGE);
    pd=fitdist(x','Normal');
% bar(BIN_RANGE,y_data./sum(y_data)/BIN_WIDTH,'facealpha',.5,'edgecolor','none')
bar(BIN_RANGE,y_data./sum(y_data)/BIN_WIDTH)

hold on
xlim([-70 -30])
x_values = -90:0.01:BIN_MAX;
hold on
y = pdf(pd,x_values);
plot(x_values,y,'LineWidth',4,'Color','r')
axis square; box off;




figure(100)
subSamp=1000;
plot(FrE(1:subSamp),CVEvec(1:subSamp),'o')
box off
ylim([0 3])
xlim([1 40])
hold all
axis square
plot(FrI(1:subSamp),CVIvec(1:subSamp),'or')
box off
ylim([0 3])
xlim([1 40])
hold all
axis square
%%

 BIN_WIDTH = 0.2;
BIN_MAX = 3;
BIN_RANGE = -4:BIN_WIDTH:BIN_MAX;
x=FrE;
mean(x)
std(x)
x(x<0.05)=[];
% x(x<1)=[];
% parmhat = lognfit(x);
pd=fitdist(log10(x),'Normal');
% Plot comparison of the histogram of the data, and the fit
figure
hold on
% Empirical distribution
y_data=hist(log10(x),BIN_RANGE);
% bar(BIN_RANGE,y_data./sum(y_data)/BIN_WIDTH,'b','facecolor','k','facealpha',.5,'edgecolor','none')
bar(BIN_RANGE,y_data./sum(y_data)/BIN_WIDTH,'b','facecolor','k','edgecolor','none')


hold on
x_values = BIN_RANGE;
y = pdf(pd,x_values);
plot(x_values,y,'LineWidth',4)
xlim([-3 3])
%
x=FrI;
mean(x)
std(x)
x(x<0.01)=[];

% x(x<1)=[];
pd=fitdist(log10(x),'Normal');
% Plot comparison of the histogram of the data, and the fit
hold on
% Empirical distribution
y_data=hist(log10(x),BIN_RANGE);
% bar(BIN_RANGE,y_data./sum(y_data)/BIN_WIDTH,'r','facecolor','k','facealpha',.5,'edgecolor','none')
bar(BIN_RANGE,y_data./sum(y_data)/BIN_WIDTH,'r','facecolor','k','edgecolor','none')

hold on
x_values = BIN_RANGE;
y = pdf(pd,x_values);
plot(x_values,y,'r','LineWidth',4)
xlim([-3 3])
%%  
figure
plot(results.V1E(1:499),FrEfromV*100,'o') 
xlim([-70 -30])
ylim([0 40])
axis square; box off
