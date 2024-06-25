% z=readtable('NMO_161366-voltageDistributionPreCosyne.csv'); % old one- on eLife folder
% z= readtable('/Users/darshanr/Dropbox (HHMI)/VoltageWholeCellProject-2018/MultiCompProject/all_chart_output_v3.csv');

 z= readtable('/Users/darshan/Dropbox/Work/VoltageWholeCellProject-2018/MultiCompProject/all_chart_output_v3_v1.csv');
%z= readtable('/Users/darshan/Dropbox/Work/VoltageWholeCellProject-2018/MultiCompProject/NMO_161366_different_location.csv');



%% Plot on-soma Vs distributed- summary
Morphology='NMO_161366'; % L5 mouse
% Morphology='NMO_130658'; % L4 FS mouse
% Morphology='NMO_02484'; % L4 Spiny Stellate mouse

% Morphology='JM072303'; % L4 Spiny Stellate mouse


for kk=1:2
    
    if kk==1
        min_dist=-0.5; max_dist=0.001; condValue=2;
    else
        min_dist=0; max_dist=9e10; condValue=2;
    end
    
    z_morph=z;
    zdis=z_morph(strcmp(z_morph.model_name,Morphology),:);
    zdis=zdis((zdis.min_dist==min_dist)&(zdis.max_dist<max_dist),:);
    
    cond_list=unique([zdis.cond_m]);
    cond_m=cond_list(condValue) % The first, which is 0.5nanoSimense
    zdis_cond=zdis(zdis.cond_m==cond_m,:);
   % MODEL NAME ADD !!!!! NMO_161366
    size(zdis_cond)
    syn_n_list=unique([zdis_cond.syn_n]);
    
    figure
    % % % syn_n_list=[204 273 409 546 819 1023 2047];
    % syn_n_list=[ 409 1023 2047 4095];
    
    % syn_n_list=[ 468 937 2344 4688];
    sig_vec=[];sdv_vec=[];totSim=[];
    for ii=1:length(syn_n_list)
        %     ii=5
        v=zdis_cond.v_mean(zdis_cond.syn_n==syn_n_list(ii));
        sig=zdis_cond.v_std(zdis_cond.syn_n==syn_n_list(ii));
        plot(v,sig,'o')
        hold on
        sig_vec(ii)=mean(sig);sdv_vec(ii)=std(v);totSim(ii)=length(v);
    end
    figure(2)
    subplot 211
    semilogx(syn_n_list,sdv_vec,'o-')
    hold on;axis square;box off;
    xlim([450 10000])
    ylabel('\Delta_v-  heteroginity voltage')
    subplot 212
    semilogx(syn_n_list,sig_vec,'o-')
    hold on;axis square;box off;
  
    
    xlim([450 10000]);ylim([0 6])
    xlabel('num of conn.')
    ylabel('std of voltage')
    % legend('on-soma,g_s=0.01','uniform-dist,g_s=0.0025')
    legend('on-soma,g_s=0.0025','uniform-dist,g_s=0.0025')
end
%% Plot on-soma Vs distributed- voltage distributions
lineStyles = linspecer(100,'sequential');
c=colormap(lineStyles);
BIN_WIDTH = 1;
BIN_MAX = -30;
BIN_RANGE = -80:BIN_WIDTH:BIN_MAX;
x_values = -90:0.01:BIN_MAX;
c_vec=[1 100];
for kk=1:2
    if kk==1
%         min_dist=-0.5; max_dist=0.001; condValue=1; numSyn=468;numSynIdx=4;

% %             min_dist=-0.5; max_dist=0.001; condValue=1; numSyn=937;numSynIdx=5;
%             min_dist=-0.5; max_dist=0.001; condValue=1; numSyn=4688;numSynIdx=7;
%             min_dist=-0.5; max_dist=0.001; condValue=2; numSyn=2344;numSynIdx=6;

            
                        min_dist=-0.5; max_dist=0.001; condValue=2; numSyn=2344;
%                       min_dist=-0.5; max_dist=0.001; condValue=2; numSyn=7813;


else
%         min_dist=0; max_dist=9e10; condValue=1; numSyn=409;numSynIdx=5;

%             min_dist=0; max_dist=9e10; condValue=1; numSyn=1023;numSynIdx=8;
%             min_dist=0; max_dist=9e10; condValue=1; numSyn=2047;numSynIdx=10;

%             min_dist=0; max_dist=9e10; condValue=1; numSyn=4095;numSynIdx=11;
            min_dist=0; max_dist=9e10; condValue=2; numSyn=2047;
%                min_dist=0; max_dist=9e10; condValue=2; numSyn=8190;

    end
    
%     zdis=z;
    z_morph=z;
    zdis=z_morph(strcmp(z_morph.model_name,Morphology),:);
    
    zdis=zdis((zdis.min_dist==min_dist)&(zdis.max_dist<max_dist),:);
    
    cond_list=unique([zdis.cond_m]);
    cond_m=cond_list(condValue);
    zdis_cond=zdis(zdis.cond_m==cond_m,:);
    size(zdis_cond)
    syn_n_list=unique([zdis_cond.syn_n])
    numSynIdx=find(syn_n_list==numSyn);
    
    figure(10)
    v=zdis_cond.v_mean(zdis_cond.syn_n==syn_n_list(numSynIdx));
    length(v)
    x=v;
    pd=fitdist(x,'Normal');
    % Plot comparison of the histogram of the data, and the fit
    % figure
    
    % Empirical distribution
    y_data=hist(x,BIN_RANGE);
    % bar(BIN_RANGE,y_data./sum(y_data)/BIN_WIDTH,'facealpha',.5,'edgecolor','none')
    bar(BIN_RANGE,y_data./sum(y_data)/BIN_WIDTH,'facealpha',.5,'edgecolor','none','facecolor',c(c_vec(kk),:))
        hold on
    % x_values = BIN_RANGE;
    y = pdf(pd,x_values);
    plot(x_values,y,'LineWidth',4,'Color',c(c_vec(kk),:))
    % bar(BIN_RANGE,y_data./sum(y_data)/BIN_WIDTH,'b','facecolor','k','edgecolor','none')
    hold on
    pause
    xlim([-70 -30])
    axis square; box off;
    hold on
       
%     legend('on-soma,g_s=0.0025,#syn=937','uniform-dist,g_s=0.0025,#syn=1023','uniform-dist,g_s=0.0025,#syn=1023')
    legend('on-soma,g_s=0.0005,#syn=2344','uniform-dist,g_s=0.0005,#syn=2047','uniform-dist,g_s=0.0025,#syn=1023')

end




