%% Make the LIF analytical figure
% 1. Generate CV, DV Vs sigma_V


% g=0.05;C=1;tau=C/g;
g_vec=[0.1 0.2 0.3]
N=1000;
% g_vec=[0.1]
for gg=1:length(g_vec)
    g=g_vec(gg);
    ;C=1;tau=C/g;
    % g=.1 ;C=1;tau=C/g;
    
    baseLine=-40;
    % % % % % % % % % % % % % % % % % % % % % % % % % sigTheta=2;
    sigTheta=0;
    
    % Vrest=baseLine-15;Theta=baseLine-0;Vr=baseLine-10;
    Vrest=baseLine-15;Theta=baseLine-0;Vr=baseLine-12;
    
    muTheta=Theta;
    
    % Vrest=-55;Theta=-40;Vr=-55;
    
    %  Vrest=-50;Theta=-37;Vr=-48; sigTheta=3;
    % Vrest=15;Theta=20;Vr=10;
    
    binSize=N/10;
    NetParm.g=g;NetParm.C=C;
    NetParm.Vrest=Vrest;NetParm.Vr=Vr;
    NetParm.sigTheta=sigTheta;
    NetParm.muTheta=muTheta;
    NetParm.N=N;
    
    
    % Network statistics
    NetParm.Theta=Theta;
    
    for it=1
        it
        t_0=0;dt=0.05;t_end=5000;
        
        switch it
            case 6
                barNu=4*(tau/1000);   barNu2=(8*tau/1000)^2;
                
            case 5
                barNu=0.04*(tau/1000);   barNu2=(0.08*tau/1000)^2;
            case 5
                barNu=0.07*(tau/1000);   barNu2=(0.1*tau/1000)^2;
                % %             t_end=10000;
            case 3
                % %             barNu=1*(tau/1000);   barNu2=(1*tau/1000)^2;
                barNu=3.6*(tau/1000);   barNu2=(4.3*tau/1000)^2;
                %             barNu=0.036;   barNu2=(0.043)^2;
                % barNu=0.032;   barNu2=(0.094)^2;
                
            case 4
                barNu=3.2*(tau/100);   barNu2=(9.4*tau/100)^2;
            case 1
                barNu=6*(tau/1000);   barNu2=(8*tau/1000)^2;
            case 2
                barNu=9.3*(tau/1000);   barNu2=(8.5*tau/1000)^2; % Inh of S1
                
        end
        
        StatVariables.barNu=barNu;StatVariables.barNu2=barNu2;
        
        SimParam.dt=dt;SimParam.t_end=t_end;SimParam.t_0=t_0;
        sig_vec=sqrt((0.5:0.5:7).^2./tau*2); % High
        col=hsv(length(sig_vec));    LW=2;
        
        for mm=1:length(sig_vec)
            mm
            sig=sig_vec(mm);
            if mm==1
                mu_out=-.6; del_out=.3;
            else
                mu_out=mu_vec(mm-1);del_out=del_vec(mm-1);
            end
            NetParm.Vr=Vr;
            [mu_out, del_out]=SolveLIFGivenSigma(NetParm,StatVariables,sig,mu_out, del_out)
            
            % For OU we should get sig^2/(2*tau*g^2).
            mu_vec(gg,mm)=mu_out;
            del_vec(gg,mm)=del_out;
            
            % Analytic\numeric
            [yN, sigv2N,FrN,FrPN]=VI_Curve_LIFver3_prim(mu_out,sig,del_out,N,NetParm);
            
            
            FrP_vec(gg,mm)=mean(FrPN); Fr_vec(gg,mm)=mean(FrN);sdFr_vec(gg,mm)=std(FrN);
            V_vec(gg,mm)=mean(yN);DelatV_vec(gg,mm)=std(yN);skV_vec(gg,mm)=skewness(yN);
            sigV_vec(gg,mm)=mean(sqrt(sigv2N));sdSDV_vec(gg,mm)=std(sqrt(sigv2N));
            VarV_vec(gg,mm)=mean((sigv2N));sdVarV_vec(gg,mm)=std((sigv2N));
            
            
            % Analytic
            [barNu(gg,mm), barNu2(gg,mm), barV(gg,mm),barV2(gg,mm), barSigV(gg,mm), barSigV2(gg,mm)]=TestSolDistLIF(mu_out, sig, del_out,NetParm);
            
            % Simulations
            I=mu_out+del_out.*randn(NetParm.N,1);
            NetParm.Vr=Vr*ones(NetParm.N,1);
            [Tout,Yout,sp]=Test_LIF_DiffAprrox_Ver2(I,sig,NetParm,SimParam);
            
            figure(1)
            tTot=size(Yout,2);
            Yout(:,1:round(tTot/20))=[];
            Tout(1:round(tTot/20))=[];
            plot(Tout,Yout(1:10,:))
            if isempty(sp)
                FrSim=0;  CV=0;
            else
                
                SpikeToTakeOff=find(sp(:,1)>Tout(round(tTot/20)),1);
                sp(1:SpikeToTakeOff,:)=[];
                [~,idx] = sort(sp(:,2)); % sort just the first column
                z = sp(idx,:);   % sort the whole matrix using the sort indices
                
                % z=sort(sp,1);
                for ii=1:NetParm.N
                    indSpike=find(z(:,2)==ii);
                    ISI=(diff(z(indSpike,1)));
                    CV(ii)=std(ISI)./mean(ISI);
                end
                
                numSpikes=histc(sp(:,2),1:NetParm.N);
                FrSim=numSpikes./(Tout(end)-Tout(round(tTot/20))).*1000/tau;
            end
            
            sim_mv=mean(Yout,2);sim_sd=std(Yout,[],2);
            
            simFr_vec(gg,mm)=mean(FrSim);simsdFr_vec(gg,mm)=std(FrSim);
            simV_vec(gg,mm)=mean(sim_mv);simsDelatV_vec(gg,mm)=std(sim_mv);
            simSigV_vec(gg,mm)=mean(sim_sd);
            
            simSKV_vec(gg,mm)=mean(skewness(Yout,[],2));
            CVplot=CV; CVplot(isnan(CVplot))=[];
            mV=mean(CVplot);SDV=std(CVplot);
            simCV_vec(gg,mm)=mV;simSdCV_vec(gg,mm)=SDV;
            % End loop
        end
    end
end
%%
tau_vec=C./g_vec;

figure(10)

% Simulation
subplot 511
plot(simSigV_vec',simsDelatV_vec','-o')
hold all;axis square;box off
title('Simulation')
subplot 512
plot(simSigV_vec',simCV_vec','-o')
hold all;axis square;box off
subplot 513
plot(sigV_vec',simV_vec','-o')
hold all;axis square;box off
subplot 514
plot(sigV_vec',simFr_vec','-o')
hold all;axis square;box off
subplot 515
plot(sigV_vec',simSKV_vec','-o')
hold all;axis square;box off

% Analytic\numeric
figure(11)
subplot 411
plot(sigV_vec',DelatV_vec','-o')
hold all;axis square;box off
title('Analytic\numeric')
subplot 412
plot(sigV_vec',simCV_vec','-o')
hold all;axis square;box off
subplot 413
plot(sigV_vec',V_vec','-o')
hold all;axis square;box off
subplot 414
plot(sigV_vec',Fr_vec'*1000./tau_vec,'-o')
hold all;axis square;box off

% Analytic

% barNu(gg,mm), barNu2(gg,mm), barV(gg,mm),barV2(gg,mm), barSigV(gg,mm), barSigV2(gg,mm)

figure(111)
subplot 411
plot(sqrt(barSigV'),sqrt(barV2'),'-o')
hold all;axis square;box off
title('Analytic')
subplot 412
plot(sqrt(barSigV'),simCV_vec','-o')
hold all;axis square;box off
subplot 413
plot(sqrt(barSigV'),barV','-o')
hold all;axis square;box off
subplot 414
plot(sqrt(barSigV'),barNu'*1000./tau_vec,'-o')
hold all;axis square;box off

figure(12)
subplot 221
plot(V_vec,simV_vec,'-o')
hold all
plot(V_vec,V_vec,'-k')
subplot 222
plot(DelatV_vec,simsDelatV_vec,'-o')
hold all
plot(DelatV_vec,DelatV_vec,'-k')
subplot 223
plot(sigV_vec,simSigV_vec,'-o')
hold all
plot(sigV_vec,sigV_vec,'-k')
subplot 224
plot(Fr_vec*1000/tau,simFr_vec,'-o')
%     plot(Fr_vec*100,simFr_vec,'-o')
hold all
plot(Fr_vec,Fr_vec,'-k')
%% 2. Generate sigma_V Vs mean V for low and high noise

% sig_vec_all=sqrt((0.5:0.5:7).^2./tau*2); % High
% sig_vec=sig_vec_all(3:2:9);

sig_vec_all=sqrt((0.5:0.05:7).^2./tau*2); % High
sig_vec=sig_vec_all(10:2:90);

        
g_vec=[0.1 0.2 0.3]
N=10000;
% g_vec=[0.1]
gg=1
barNu=6*(tau/1000);   barNu2=(8*tau/1000)^2;
g=g_vec(gg);
C=1;tau=C/g;
% g=.1 ;C=1;tau=C/g;

baseLine=-40;
% % % % % % % % % % % % % % % % % % % % % % % % % sigTheta=2;
sigTheta=0;

% Vrest=baseLine-15;Theta=baseLine-0;Vr=baseLine-10;
Vrest=baseLine-15;Theta=baseLine-0;Vr=baseLine-12;

muTheta=Theta;

% Vrest=-55;Theta=-40;Vr=-55;

%  Vrest=-50;Theta=-37;Vr=-48; sigTheta=3;
% Vrest=15;Theta=20;Vr=10;

binSize=N/10;
NetParm.g=g;NetParm.C=C;
NetParm.Vrest=Vrest;NetParm.Vr=Vr;
NetParm.sigTheta=sigTheta;
NetParm.muTheta=muTheta;
NetParm.N=N;


% Network statistics
NetParm.Theta=Theta;


StatVariables.barNu=barNu;StatVariables.barNu2=barNu2;

SimParam.dt=dt;SimParam.t_end=t_end;SimParam.t_0=t_0;
col=hsv(length(sig_vec));    LW=2;

for mm=1:length(sig_vec)
    mm
    sig=sig_vec(mm);
    if mm==1
        mu_out=-.6; del_out=.3;
    else
        mu_out=mu_vec(mm-1);del_out=del_vec(mm-1);
    end
    NetParm.Vr=Vr;
    [mu_out, del_out]=SolveLIFGivenSigma(NetParm,StatVariables,sig,mu_out, del_out)
    
    % For OU we should get sig^2/(2*tau*g^2).
    mu_vec(mm)=mu_out;
    del_vec(mm)=del_out;
    
    % Analytic\numeric
    [yN, sigv2N,FrN,FrPN]=VI_Curve_LIFver3_prim(mu_out,sig,del_out,N,NetParm);
    
    figure(5)
    plot(yN,sqrt(sigv2N),'o'); hold on
    
    cor_vec(mm)=corr(yN,sqrt(sigv2N));
    
    box off; axis square; xlim([-80 -30]); ylim ([0 6])
    sigV(mm)=mean(sqrt(sigv2N));
end
legend(sprintf('sigma=%.1f;  sigma_v=%.1f',sig_vec(1),sigV(1)),  sprintf('sigma=%.1f;  sigma_v=%.1f',sig_vec(2),sigV(2)), sprintf('sigma=%.1f;  sigma_v=%.1f',sig_vec(3),sigV(3)), sprintf('sigma=%.1f;  sigma_v=%.1f',sig_vec(4),sigV(4)))
figure(6)
plot(sigV,cor_vec,'-')
xlim([2 4]); ylim([-1 1]); box off; axis square
