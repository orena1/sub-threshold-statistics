%%
d=1.5*1*10^(-6);%(in microm);
g_list =[0.1000 0.5000 1.0000 2.5000 5.0000];
    
CellNameCell{1}='L5MouseR1';
CellNameCell{2}='BasketR1';
CellNameCell{3}='SpinyR1';


for  ii=1:1
    CellName=CellNameCell{ii};
for gg=2

% Layer 5 mouse
if strcmp(CellName,'L5MouseR1')
z= readtable(sprintf('/Users/darshan/Dropbox/Work/VoltageWholeCellProject-2018/MultiCompProject/rhos/Different_g_syn/NMO_161366/Rm3500_Rin57_norm_%d.csv',gg));
end

% Layer 4 Spiny
if strcmp(CellName,'SpinyR1')
z= readtable(sprintf('/Users/darshanr/Dropbox (HHMI)/VoltageWholeCellProject-2018/MultiCompProject/rhos/Different_g_syn/JM072303/Rm1000_Rin45_norm_%d.csv',gg));
end


if strcmp(CellName,'BasketR1')
z= readtable(sprintf('/Users/darshanr/Dropbox (HHMI)/VoltageWholeCellProject-2018/MultiCompProject/rhos/Different_g_syn/NMO_130658/Rm3000_Rin36_norm_%d.csv',gg));
end


if ii==1
Rm=3500; Rin=57;
end
if ii==2
    Rm=7000;Rin=94;
end
if ii==3
    Rm=14000; Rin=162;
end

lambda(ii)=(d*Rm/(4*Rin))^0.5 *10e4 %results: 612 microM





[rhoE,alphaE]=ChooseAlphaRhoFromDataEI(z,1000);


x=z.seg_dist(2:end);
alpha=z.alpha(2:end);
rho=z.rho(2:end);
figure(1)
subplot 311
plot(x,rho,'.')
hold on;box off; ylim([0 1]);
ylabel('rho')
subplot 312
plot(x,alpha,'.')
hold on;box off; ylim([0 1]);
ylabel('alpha')
subplot 313
plot(x,alpha.*rho,'.')
hold on;box off; ylim([0 1]);
ylabel('Gamma')

figure(111)
semilogy(x,rho,'.')
hold on;box off; ylim([0 1]);
ylabel('rho')

% figure(112)
% semilogy(x,alpha,'.')
% hold on;box off; ylim([0 1]);
% ylabel('alpha')

figure(100*ii+12)
plot(x,log10(alpha),'.')
hold on;box off; 
ylabel('alpha')
y=alpha;
[RsqNW,yCalc]= MyLinearReg(x,log10(y));
plot(x,yCalc,'r')
[RsqNW,yCalc,expAlpha(ii,:)]= MyLinearReg(x,log(y));
ylim([-3 0]);


figure(100*ii+13)
plot(x,log10(alpha.*rho),'.')
hold on;box off; 
% ylim([0 1]);
ylabel('Gamma')
y=alpha.*rho;
[RsqNW,yCalc]= MyLinearReg(x,log10(y));
plot(x,yCalc,'r')
ylim([-3 0]);
[RsqNW,yCalc,expGamma(ii,:)]= MyLinearReg(x,log(y));



% 
% figure(2)
% subplot 211
% plot(x./lambda(ii),rho,'.')
% hold on;box off;axis square; ylim([0 1]);
% ylabel('rho')
% subplot 212
% plot(x./lambda(ii),alpha,'.')
% hold on;box off;axis square; ylim([0 1]);
% ylabel('alpha')

figure(3)
plot(alphaE,rhoE.*alphaE,'.')
hold on;box off;axis square; xlim([0 1]);  ylim([0 1]);

xlabel('Somatic current change (\alpha)')
ylabel('Visability of conductance (\alpha\rho)')
m_rho(ii,gg)=mean(rhoE);
m_alpha(ii,gg)=mean(alphaE);
m_mult(ii,gg)=mean(rhoE.*alphaE);


% Distributions

figure(100+ii)
histogram2(rhoE,alphaE,0:0.05:1,0:0.05:1,'DisplayStyle','tile','ShowEmptyBins','on');
title('rho')
xlim([0 1]);ylim([0 1])
axis square

figure(200+ii)
hist(alphaE,0:0.05:1)
title('alpha')
xlim([0 1]);
box off

figure(300+ii)
hist(rhoE,0:0.05:1)
title('rho')
xlim([0 1]);
box off


end

end
figure(1)
legend(['L5MouseR1, \langle\rho\rangle=' num2str(m_rho(1)) '  \langle\alpha\rangle=' num2str(m_alpha(1))],['BasketR1, \langle\rho\rangle=' num2str(m_rho(2)) '  \langle\alpha\rangle=' num2str(m_alpha(2))],['SpinyR1, \langle\rho\rangle=' num2str(m_rho(3)) '  \langle\alpha\rangle=' num2str(m_alpha(3))])
figure(3)

legend(['L5MouseR1, \langle\alpha\rho\rangle=' num2str(m_mult(1)) '  \langle\alpha\rangle=' num2str(m_alpha(1))],['BasketR1, \langle\alpha\rho\rangle=' num2str(m_mult(2)) '  \langle\alpha\rangle=' num2str(m_alpha(2))],['SpinyR1, \langle\alpha\rho\rangle=' num2str(m_mult(3)) '  \langle\alpha\rangle=' num2str(m_alpha(3))])
alphavec=0:0.01:1;
plot(alphavec,alphavec,'--k')
plot(alphavec,alphavec.^2,'-r')
%%
close all
for ii=1
figure(100+ii)
loglog(g_list,m_mult(ii,:),'-o')
hold all
plot(g_list,m_alpha(ii,:),'-o')
plot(g_list,m_rho(ii,:),'-o')
loglog(g_list,1-0.1*g_list,'-ok')
 xlabel('g_syn')
 legend('\langle \Gamma \rangle=\langle\alpha\rho\rangle','\langle\alpha\rangle','\langle\rho\rangle')
end


