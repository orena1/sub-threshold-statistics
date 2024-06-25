% Load data
 
Path='.../ALM/SiliconProbeData/FixedDelayTask';
cd (Path)



RFlag=1;

fileListUnSorted = getAllFiles(Path);
fileListUnSorted(end)=[];
fileListUnSorted(1)=[];
fileListUnSorted(end)=[];

fileList=natsortfiles(fileListUnSorted);
numFiles=length(fileList);
TotalUnitsAllFiles=0;
uu=0;
PSTHD=[];
for ff=1:length(fileList)
%     for ff=1

    ff
    fileName = fileList{ff};
    % unitNum  = 25;
    preCueDur  = 4.2; % in sec: plot how many seconds before go cue onset.
    postCueDur = 1.4; % in sec: plot how many seconds after go cue onset.
    
    
    
    % laod file & extract info required for plotting
    load(fileName)
    TotalUnits=length(unit);
    TotalUnitsAllFiles=TotalUnitsAllFiles+TotalUnits;
    
    % extract info
    % for unitNum=3
            trialType      = unit(1).Behavior.Trial_types_of_response_vector; % trial type based on outcome
        NumCorrectTrials(ff)=sum(trialType==1)+sum(trialType==2) 
        NumIncorrectTrials(ff)=sum(trialType==3)+sum(trialType==4);
        NumAllTrials(ff)=length(trialType);
        
    for unitNum=1:TotalUnits
        uu=uu+1;
        SpikeWidth=unit(unitNum).SpikeWidth;
        Depth=unit(unitNum).Depth;
        
        trialType      = unit(unitNum).Behavior.Trial_types_of_response_vector; % trial type based on outcome
 
        % 1: correct lick R trial,  2: correct lick L trial,  3: incorrect lick R trial,  4: incorrect lick L trial,
        % 5: correct lick R trial with early lick,   6: correct lick L trial with early lick,
        % 7: incorrect lick R trial with early lick, 8: incorrect lick L trial with early lick
        % 9: lick R trial but no response,          10: lick L trial but no response,
        % 11: others(unidentified)
        
        photoStimTrial = unit(unitNum).Behavior.stim_trial_vector;  % photo stim trial type
        % 0: no stim,   1: 0.05mW bilateral stim during early dealy,
        % 2: 0.1mW bilateral stim during early dealy,
        % 3: 0.2mW bilateral stim during early dealy,
        % 4: 0.3mW bilateral stim during early dealy,
        % don't analyze 0.05mW it was too weak to have any behavioral effect
        
        
        % extract timing info
        T.tAxisPSTH   = -preCueDur-0.1:0.001:postCueDur+0.1;      % T axis for PSTH, shift by 0.1 to remove smoothing artifact at the edge
        T.sampleOnset = [unit(unitNum).Behavior.Sample_start];  % smaple epoch onset
        T.delayOnset  = [unit(unitNum).Behavior.Delay_start];     % delay epoch onset
        T.cueOnset    = [unit(unitNum).Behavior.Cue_start];       % go cue  onset
        
        % calculate mean onset time of each epoch for plotting
        % exclude early lick and non-response trial by "trialType<5"
        T.meanSampleOnset = mean(T.sampleOnset(trialType==1));
        T.meanDelayOnset  = mean(T.delayOnset(trialType==1));
        T.meanCueOnset    = mean(T.cueOnset(trialType==1));
        
        spikeTimes       = unit(unitNum).SpikeTimes; % time bin of spike
        trialIdxOfSpikes = unit(unitNum).Trial_idx_of_spike; % time bin of spike
        
        trialRange            = unit(unitNum).Trial_info.Trial_range_to_analyze; % range of trials to analyze
        trialTypeInRange      = trialType(trialRange(1) : trialRange(2));        % trial type in range
        photoStimTrialInRange = photoStimTrial(trialRange(1) : trialRange(2));   % stim trial type in range
        
        
        selectedTrials=(photoStimTrialInRange==0).*(trialTypeInRange==RFlag);
        FrD=[];      FrS=[];CV2S=[];CVS=[];
        ISID=[];ISIS=[];

        %         for tr = trialRange(1) : trialRange(2)
        trVec=trialRange(1) : trialRange(2);
        
        trInd=0;
        for tr = trVec(selectedTrials==1)
            trInd=trInd+1;
            
            % extarct spikes of each trial
            %             SpikesTmp   = spikeTimes(trialIdxOfSpikes == tr) - T.cueOnset(tr); % align to go cue
            sampleInd=(spikeTimes(trialIdxOfSpikes == tr)<T.delayOnset(tr)).*(spikeTimes(trialIdxOfSpikes == tr)>T.sampleOnset(tr));
            delayInd=(spikeTimes(trialIdxOfSpikes == tr)<T.cueOnset(tr)).*(spikeTimes(trialIdxOfSpikes == tr)>T.delayOnset(tr));
            cueInd=(spikeTimes(trialIdxOfSpikes == tr)>T.cueOnset(tr));
            
            sampleTot=T.delayOnset(tr)-T.sampleOnset(tr);
            delayTot=T.cueOnset(tr)-T.delayOnset(tr);
            cueTot=postCueDur;
            
%             SpikesTmpD   = spikeTimes(trialIdxOfSpikes == tr) - T.delayOnset(tr); % align to go cue

            SpikesTmp   = spikeTimes(trialIdxOfSpikes == tr) ; 
            SpikesD=SpikesTmp(delayInd==1);
%             PSTHD{uu}=[PSTHD{uu} SpikesD tr*length(SpikesD)];
            SpikesS=SpikesTmp(sampleInd==1);
            SpikesR=SpikesTmp(cueInd==1);
            
            if length(SpikesTmp)<2
                FrD(trInd)=0;
                CVD(trInd)=0;
                CV2D(trInd)=0;
                
                FrS(trInd)=0;
                CVS(trInd)=0;
                CV2S(trInd)=0;
                
                FrR(trInd)=0;
            else
                FrS(trInd)=length(SpikesS)./(sampleTot);
                FrD(trInd)=length(SpikesD)./(delayTot);
                  FrR(trInd)=length(SpikesR)./(cueTot);
                
                % Calc CV
                if sum(FrD>60)
                    unitNum
                    trInd
                end
                ISID=[ISID; diff(SpikesD)];
                                
                ISIS=[ISIS; diff(SpikesS)];
                        
            end
    
        end   
        Session(ff).FrDUnit(unitNum)=mean(FrD);
        Session(ff).FrRUnit(unitNum)=mean(FrR);
        Session(ff).CVDUnit(unitNum)=std(ISID)/mean(ISID);
        Session(ff).CVDUnitTrials(unitNum)=length(ISID);
         Session(ff).CVSUnitTrials(unitNum)=length(ISIS);
        Session(ff).ISID{unitNum}=ISID;
        Session(ff).ISIS{unitNum}=ISIS;

        Session(ff).CV2DUnit(unitNum)=2.*mean(abs(ISID(2:end)-ISID(1:end-1))./(ISID(2:end)+ISID(1:end-1)));
        Session(ff).FrSUnit(unitNum)=mean(FrS);
        Session(ff).CVSUnit(unitNum)=std(ISIS)/mean(ISIS);
        Session(ff).CV2SUnit(unitNum)=2.*mean(abs(ISIS(2:end)-ISIS(1:end-1))./(ISIS(2:end)+ISIS(1:end-1)));
        Session(ff).SpikeWidth(unitNum)=SpikeWidth;
        Session(ff).Depth(unitNum)=Depth;
        
        
        
    end
    

    
    FS_vec=Session(ff).SpikeWidth<0.35;
    Exc_vec=Session(ff).SpikeWidth>0.5;
    
%     Layer5=(Session(ff).Depth>400) .* (Session(ff).Depth<1100);
    Layer5=(Session(ff).Depth>0) .* (Session(ff).Depth<1300);
    Exc_vecLayer5=logical(Layer5.*Exc_vec);
    figure(1000)
    plot(Session(ff).FrDUnit(Exc_vecLayer5), Session(ff).Depth(Exc_vecLayer5),'ob')
     plot(Session(ff).FrDUnit(FS_vec), Session(ff).Depth(FS_vec),'og')
    hold on
    box off
    
    figure(10)
    plot(Session(ff).FrSUnit(Exc_vec),Session(ff).CVSUnit(Exc_vec),'or')
    hold on
    plot(Session(ff).FrDUnit(Exc_vec),Session(ff).CVDUnit(Exc_vec),'ob')
    box off
    ylim([0 3])
    xlim([0 70])
    axis square
    
    figure(110)
    plot(Session(ff).FrSUnit(FS_vec),Session(ff).CVSUnit(FS_vec),'or')
    hold on
    plot(Session(ff).FrDUnit(FS_vec),Session(ff).CVDUnit(FS_vec),'og')
    box off
    ylim([0 3])
    xlim([0 70])
    axis square
    
    
    figure(12)
    plot(Session(ff).FrSUnit(Exc_vec),Session(ff).CV2SUnit(Exc_vec),'or')
    hold on
    plot(Session(ff).FrDUnit(Exc_vec),Session(ff).CV2DUnit(Exc_vec),'ob')
    box off
    ylim([0 3])
    xlim([0 70])
    axis square
    
    figure(112)
    plot(Session(ff).FrSUnit(FS_vec),Session(ff).CV2SUnit(FS_vec),'or')
    hold on
    plot(Session(ff).FrDUnit(FS_vec),Session(ff).CV2DUnit(FS_vec),'og')
    box off
    ylim([0 3])
    xlim([0 70])
    axis square

    Session(ff).fileName=fileName;

end

%%
EtotFrR=[];EtotFrS=[];EtotFrD=[];EtotCVS=[];EtotCVD=[];EtotCV2S=[];EtotCV2D=[];
ItotFrS=[];ItotFrD=[];ItotCVS=[];ItotCVD=[];ItotCV2S=[];ItotCV2D=[];
EtotFrDL5=[];EtotCVDL5=[];DepthAll=[];EtotCVDTrials=[];
Exc_vec_List=[];
FS_vec_List=[];
for ff=1:numFiles
    FS_vec=Session(ff).SpikeWidth<0.35;
    Exc_vec=Session(ff).SpikeWidth>0.5;
    Layer5=(Session(ff).Depth>400) .* (Session(ff).Depth<800);
    Exc_vecLayer5=logical(Layer5.*Exc_vec);
    
    DepthAll=[DepthAll Session(ff).Depth(Exc_vec)];
    
    EtotFrDL5=[EtotFrDL5 Session(ff).FrDUnit(Exc_vecLayer5)];
    EtotCVDL5=[EtotCVDL5 Session(ff).CVDUnit(Exc_vecLayer5)];
Exc_vec_List=[Exc_vec_List Exc_vec];
FS_vec_List=[FS_vec_List FS_vec];
        
    EtotFrS=[EtotFrS Session(ff).FrSUnit(Exc_vec)];
    EtotFrD=[EtotFrD Session(ff).FrDUnit(Exc_vec)];
        EtotFrR=[EtotFrR Session(ff).FrRUnit(Exc_vec)];
    EtotCVS=[EtotCVS Session(ff).CVSUnit(Exc_vec)];
    EtotCVD=[EtotCVD Session(ff).CVDUnit(Exc_vec)];

     EtotCVDTrials=[EtotCVDTrials Session(ff).CVDUnitTrials(Exc_vec)];

    EtotCV2S=[EtotCV2S Session(ff).CV2SUnit(Exc_vec)];
    EtotCV2D=[EtotCV2D Session(ff).CV2DUnit(Exc_vec)];
    
    ItotFrS=[ItotFrS Session(ff).FrSUnit(FS_vec)];
    ItotFrD=[ItotFrD Session(ff).FrDUnit(FS_vec)];
    ItotCVS=[ItotCVS Session(ff).CVSUnit(FS_vec)];
    ItotCVD=[ItotCVD Session(ff).CVDUnit(FS_vec)];
    ItotCV2S=[ItotCV2S Session(ff).CV2SUnit(FS_vec)];
    ItotCV2D=[ItotCV2D Session(ff).CV2DUnit(FS_vec)];
end
figure(100)
hist(EtotFrD,25)
box off
axis square
xlim([0 70])

CVD=EtotCVD;
CVD(isnan(EtotCVD))=[];
CVD(CVD==0)=[];
figure(200)
hist(CVD,25)
box off
xlim([0 3])
axis square


figure(600)
% subplot 211
plot(EtotFrD,EtotCVD,'o')
hold on
plot(ItotFrD,ItotCVD,'o')
box off
axis square
ylim([0 3])
%%
figure(222)
plot(EtotCVS,EtotCVD,'o')
hold on
plot(ItotCVS,ItotCVD,'or')
plot(0:1:3,0:1:3,'-k')
box off
axis square
ylim([0 3])
xlim([0 3])
EtotCVD(isnan(EtotCVD))=[];
EtotCVS(isnan(EtotCVS))=[];
[HSD,PSD,CI,STATS] = ttest2(EtotCVS,EtotCVD);

%% Figure 1C: Log normal (log10) E and I
% Number of data points
N = 10000;
% Specify bin locations for histogram and fit
BIN_WIDTH = 0.1;
BIN_MAX = 2;
BIN_RANGE = -2:BIN_WIDTH:BIN_MAX;
% x_values = BIN_RANGE;
% Make up some data. (You should use your real data in place of x.)
x = EtotFrD';%lognrnd(1,0.3,N,1);
% Fit the data
x(find(EtotFrD==0))=0.001;
% parmhat = lognfit(x);
pd=fitdist(log10(x),'Normal');
% Plot comparison of the histogram of the data, and the fit
figure
hold on
% Empirical distribution
y_data=hist(log10(x),BIN_RANGE);
y_clear=hist((x),0:0.1:60);
bar(BIN_RANGE,y_data./sum(y_data)/0.1,'b','facecolor','k','facealpha',.5,'edgecolor','none')
y = pdf(pd,BIN_RANGE);
plot(BIN_RANGE,y,'b','LineWidth',4)
% semilogx(0:0.1:60,y_clear)
% hist(log10(x)./sum(log10(x))/0.1,BIN_RANGE)

hold on
BIN_WIDTH = 0.2;
BIN_MAX = 3;
BIN_RANGE = -3:BIN_WIDTH:BIN_MAX;


x = ItotFrD';%lognrnd(1,0.3,N,1);
% Fit the data
x(find(ItotFrD==0))=0.001;
% parmhat = lognfit(x);
pd=fitdist(log10(x),'Normal');
% Plot comparison of the histogram of the data, and the fit
% figure
hold on
% Empirical distribution
y_data=hist(log10(x),BIN_RANGE);
y_clear=hist((x),0:0.1:60);
bar(BIN_RANGE,y_data./sum(y_data)/BIN_WIDTH,'r','facecolor','k','facealpha',.5,'edgecolor','none')
% semilogx(0:0.1:60,y_clear)
% hist(log10(x)./sum(log10(x))/0.1,BIN_RANGE)

hold on

x_values = BIN_RANGE;
y = pdf(pd,x_values);
plot(x_values,y,'r','LineWidth',4)
%% Figure 1F-E-G
%First load:
Path='...ALM/WholeCellData/Data';

cd (Path)
% % % % RFlag=1;
LRFlag=1; % 1=right, 2=left

skewns = @(x) (sum((x-mean(x)).^3)./length(x)) ./ (var(x,1).^1.5);
kurtss = @(x) (sum((x-mean(x)).^4)./length(x)) ./ (var(x,1).^2);
ExKurtss = @(x) (sum((x-mean(x)).^4)./length(x)) ./ (var(x,1).^2)-3;

fileListUnSorted = getAllFiles(Path);
fileListUnSorted{1}
% fileListUnSorted(1)=[];
fileList=natsortfiles(fileListUnSorted);

% % % % % % % % % % % % % fileList(1)=[];
% % % % % % % % % % % % % fileList(17)=[];
% % % % % % % % % % % % % fileList(23)=[];
Vthparams= [-55 -45 100 5];


RepCells=[7:8 12:13 29:30 33:34 38:39 47:48 50:51 52:53 54:55 56:57]; %10 cells are double
goodCellsNoRepRegular=[1:3 5:6 8 9:11 13 16:21 23:24 27 30 32 34 35:37 40:42 44 48:49 51 55 57 62 64:67 69:70 81:85 89];
% selecGoodCelssHid=selecVecHideRep(goodCellsNoRepRegular);

numFiles=size(fileList,1);
AuditoryFlag=zeros(numFiles,1);
AuditoryFlag(1:42)=1; % file 42= cell 133, afterwards it is tactile
RepCellsIndx=[7 12 29 33 38 47 50 52 54 56];
CtxDepth=[360 454 398 579 458 538 566 552 684 731 679 725 684 818 235 280 517 582 611 367 289 262 566 629 498 553 513 627 814 435 502 617 558 641 575 293 620 571 353 454 509 357 442 518 557 389 593 481 595 394 391 489 489 456 478 452 560 579 604 615 584 527 377 521 705 500 465 512 477 561 591 427 502 435 760 509 470 617 530];
CtxDepthRep=[360 454 398 579 458 538 566 566 552 684 731 679 679 725 684 818 235 280 517 582 611 367 289 262 566 629 498 553 513 513 627 814 435 435 502 617 558 641 641 575 293 620 571 353 454 509 357 357 442 518 518 557 557 389 389 593 593 481 595 394 391 489 489 456 478 452 560 579 604 615 584 527 377 521 705 500 465 512 477 561 591 427 502 435 760 509 470 617 530];

ind=0;
%
for  ff=goodCellsNoRepRegular
    tic;
    
    ind=ind+1;
    ff
    %     [data, trial, general] = JustloadNWB_Inagaki(fileList{ff});
    load(fileList{ff});
    
    % Plot
    
    % basic parameters for plotting
    %     sampleRate       = 1/mean(diff(data.tAxis )); % sample rate (per s) it should be 20kHz
    % sample rate
    sampleRate     = wholeCell.recording_data.sample_rate
    
    trialType      = wholeCell.behavioral_data.trial_type_vector; % trial type
    
    ExcludeStimTrials=strncmpi(wholeCell.behavioral_data.original_trial_name,'r_s',3)+strncmpi(wholeCell.behavioral_data.original_trial_name,'l_s',3);

    PhotoStimTrials= wholeCell.behavioral_data.AOM_on_or_off ;
    
    ExcludeTrials = (ExcludeStimTrials+PhotoStimTrials+(trialType~=LRFlag))>0;
    
    % 1: correct lick R trial,  2: correct lick L trial,  3: incorrect lick R trial,  4: incorrect lick L trial,
    % 5: correct lick R trial with early lick,   6: correct lick L trial with early lick,
    % 7: incorrect lick R trial with early lick, 8: incorrect lick L trial with early lick
    % 9: lick R trial but no response,          10: lick L trial but no response,
    % 11: others(unidentified)
    
    %     plotDuration     = 5; % 5s
    %     preDelayDuration = 2; % duration to plot before Delay onset
    preCueDur        = 3;
    postCueDur        = 1.4;
    %     numPoint         = plotDuration*sampleRate; % number of data points to plot
    %     Time.tAxis       = (1:numPoint)/sampleRate - preDelayDuration; % time axis for plotting
    
    % extract timing info
    T.tAxis       = (-preCueDur*sampleRate:postCueDur*sampleRate)/sampleRate; % T axis for Vm
    %     T.tAxisPSTH   = -preCueDur-0.1:0.001:postCueDur+0.1;                      % T axis for PSTH (add 100ms because of smoothing)
    T.sampleOnset = [wholeCell.behavioral_data.behav_timing.sampling_start];  % smaple epoch onset
    T.delayOnset  = [wholeCell.behavioral_data.behav_timing.delay_start];     % delay epoch onset
    T.cueOnset    = [wholeCell.behavioral_data.behav_timing.cue_start];       % go cue  onset
    
    
    % calculate mean onset time of each epoch for plotting
    % exclude early lick and non-response trial by "trialType<5"
    T.meanSampleOnset = mean(T.sampleOnset(trialType<5));
    T.meanDelayOnset  = mean(T.delayOnset(trialType<5));
    T.meanCueOnset    = mean(T.cueOnset(trialType<5));
    
    
    numOfTrials   = numel(wholeCell.behavioral_data.trial_type_vector); % total number of trials
    
    
    Vm          = nan(numOfTrials,numel(T.tAxis));
    Vm_wo_spike = nan(numOfTrials,numel(T.tAxis));
    spikePeaksBin = wholeCell.recording_data.spike_peak_bin; % time bin of spike
    
    
    
if AuditoryFlag(ff)
    EnResponse=find(T.tAxis==0.2);
    EnDelay=find(T.tAxis==0);
    StDelay=find(T.tAxis>-1.2,1);
    StSample=find(T.tAxis>-1.2-1.15,1);
    EnSample=StDelay;
    durSample=1.15;

    LateDelay=find(T.tAxis>-0.1,1);
else
    EnResponse=find(T.tAxis==0.2);
    EnDelay=find(T.tAxis==0);
    StDelay=find(T.tAxis>-1.2,1);
    StSample=find(T.tAxis>-1.2-1.4,1);
    EnSample=StDelay;
    durSample=1.4;
    
    LateDelay=find(T.tAxis>-0.1,1);
end
    
    Time=T.tAxis;
    Vtot=[];  SVtot=[]; DVtot=[];RVtot=[];VtotNoMeanRed=[];VtotNoMeanRedMY=[];VtotMY=[];
    SVtotMY=[]; DVtotMY=[];RVtotMY=[];
    SSK=[]; SKur=[]; Smu=[]; SSD=[];Skstat=[];
    DSK=[]; DKur=[]; Dmu=[]; DSD=[];Dkstat=[];
    RSK=[]; RKur=[]; Rmu=[]; RSD=[];
    
    SSKMY=[]; SKurMY=[]; SmuMY=[]; SSDMY=[];SkstatMY=[];
    DSKMY=[]; DKurMY=[]; DmuMY=[]; DSDMY=[];DkstatMY=[];
    RSKMY=[]; RKurMY=[]; RmuMY=[]; RSDMY=[];RkstatMY=[];
    
    SK=[]; Kur=[]; mu=[]; SD=[];kstat=[];
    files{ind}.Dmu=[];files{ind}.Smu=[];files{ind}.Rmu=[];
    files{ind}.DFr=[];files{ind}.SFr=[];files{ind}.RFr=[];
    files{ind}.RISI=[];files{ind}.DISI=[];files{ind}.SISI=[];
    ISItmp=[];
    for tr=1: numOfTrials
        % calculate the first and last bin of trial
        trialOnsetBin    = wholeCell.behavioral_data.trial_onset_bin(tr); % bin of trial onset
        cueOnset = T.cueOnset(tr); % go cue onset of each trial
        
        firstBin = trialOnsetBin + round((cueOnset-preCueDur)*sampleRate);
        lastBin  = trialOnsetBin + round((cueOnset+postCueDur)*sampleRate);
        
        
        % extarct Vm between the first and last bin
        Vm(tr,:)          = wholeCell.recording_data.Vm(firstBin:lastBin);
        Vm_wo_spike(tr,:) = wholeCell.recording_data.Vm_wo_spike(firstBin:lastBin);
    end
    
    Vm(ExcludeTrials,:)=[];
    Vm_wo_spike(ExcludeTrials,:)=[];
    numOfTrialsNotExcluded=numOfTrials-sum(ExcludeTrials)
    AllThreshVec=[];DThreshVec=[];SThreshVec=[];RThreshVec=[];
    
    TotTreshS=[];TotTreshD=[];TotTreshR=[];TotTreshAll=[];

     Trials(ind)=numOfTrialsNotExcluded;
    for tr = 1 : numOfTrialsNotExcluded

%         RVmnoap=Vm_wo_spike(tr,EnDelay:EnResponse); Change
% to Late delay in order to understand the selectivity level
 RVmnoap=Vm_wo_spike(tr,LateDelay:EnDelay);
        DVmnoap=Vm_wo_spike(tr,StDelay:EnDelay);
        SVmnoap=Vm_wo_spike(tr,StSample:EnSample);
        ALLVm=Vm(tr,:);
        
% % % % % % %         winBefore=15; USED FOR THE NEURON VERSION. NO DIFFERENCE
% % % % % % %         winAfter=15;
                winBefore=0.5;
        winAfter=10;
        
% % % % % % % % % % % % % %         RVm=Vm(tr,EnDelay:EnResponse); Change
% to Late delay in order to understand the selectivity level
        RVm=Vm(tr,LateDelay:EnDelay);

[vmtmp, Rspk ,RVthresh]=removeAPnewP_Ver2(double(RVm)',sampleRate,0.33,  Vthparams, 10000, winBefore, winAfter);
        files{ind}.RFr(tr)=sum(Rspk)/0.2; % Response Fr
        RVmnoapMY=sgolayfilt(vmtmp, 3, 21)'; % remove AP
        
        RThreshVec=[RThreshVec RVthresh{1,1}];
        RThreshVec(RThreshVec<-45)=[];% Take care of the current injections
        files{ind}.RThresholdsVec{tr}=RThreshVec;
        
        [files{ind}.RFrb(tr),RISI,files{ind}.RCV(tr),files{ind}.RCV2(tr)]=SpikeStatisticALM(sampleRate,(Rspk),1.2);
        files{ind}.RISI=[files{ind}.RISI; RISI];
        files{ind}.Rspk{tr}=Rspk;
        
        
        
        DVm=Vm(tr,StDelay:EnDelay);
        [vmtmp, Dspk ,DVthresh]=removeAPnewP_Ver2(double(DVm)',sampleRate,0.33,  Vthparams, 10000, winBefore, winAfter);
        files{ind}.DFr(tr)=sum(Dspk)/1.2; % delay Fr
        DVmnoapMY=sgolayfilt(vmtmp, 3, 21)'; % remove AP
        
%         DThreshVec=[DThreshVec DVthresh{1,1}];  % THIS SEEMS LIKE A BUG 01/09/21; Take off
        DThreshVec=DVthresh{1,1}; 

        DThreshVec(DThreshVec<-45)=[];% Take care of the current injections
        files{ind}.DThresholdsVec{tr}=DThreshVec;
        
        [files{ind}.DFrb(tr),DISI,files{ind}.DCV(tr),files{ind}.DCV2(tr)]=SpikeStatisticALM(sampleRate,(Dspk),1.2);
        files{ind}.DISI=[files{ind}.DISI; DISI];
        files{ind}.Dspk{tr}=Dspk;
        
        SVm=Vm(tr,StSample:EnSample);
        [vmtmp, Sspk,SVthresh ]=removeAPnewP_Ver2(double(SVm)',sampleRate,0.33,  Vthparams, 10000, winBefore, winAfter);
        files{ind}.SFr(tr)=sum(Sspk)/durSample; % sample Fr
        SVmnoapMY=sgolayfilt(vmtmp, 3, 21)'; % remove AP
        
%         SThreshVec=[SThreshVec SVthresh{1,1}]; % THIS SEEMS LIKE A BUG 01/09/21; Take off
SThreshVec=[ SVthresh{1,1}]; %
        SThreshVec(SThreshVec<-45)=[];% Take care of the current injections
        files{ind}.SThresholdsVec{tr}=SThreshVec;
        
        [files{ind}.SFrb(tr),SISI,files{ind}.SCV(tr),files{ind}.SCV2(tr)]=SpikeStatisticALM(sampleRate,(Sspk),1.2);
        files{ind}.SISI=[files{ind}.SISI; SISI];
        files{ind}.Sspk{tr}=Sspk;
%            ISItmp=[ISItmp;ISI];
           
        
        files{ind}.Rmu(tr)=mean(RVmnoap);
        files{ind}.Dmu(tr)=mean(DVmnoap);
        files{ind}.Smu(tr)=mean(SVmnoap);
        [vmtmp, spkMY, ALLthresh ]=removeAPnewP_Ver2(double(ALLVm)',sampleRate,0.33,  Vthparams, 10000, winBefore, winAfter);
        VmnoapMY=sgolayfilt(vmtmp, 3, 21); % remove AP
        
        [files{ind}.ALLFr(tr),ISI,files{ind}.ALLCV(tr),files{ind}.ALLCV2(tr)]=SpikeStatisticALM(sampleRate,(spkMY),5);
        files{ind}.spkMY{tr}=spkMY;
        
%         ISItmp=[ISItmp;ISI];
        
%         AllThreshVec=[AllThreshVec ALLthresh{1,1}]; % THIS SEEMS LIKE A BUG 01/09/21; Take off
        AllThreshVec=[ ALLthresh{1,1}];
        AllThreshVec(AllThreshVec<-45)=[];% Take care of the current injections
        files{ind}.AllThresholdsVec{tr}=AllThreshVec;
        
        TotTreshR=[TotTreshR RThreshVec];
        TotTreshS=[TotTreshS SThreshVec];
        TotTreshD=[TotTreshD DThreshVec];
        TotTreshAll=[TotTreshAll AllThreshVec];
        
        
        Vmnoap=Vm_wo_spike(tr,:);
        
        
        SSK(tr)=skewns(SVmnoap) ;
        SKur(tr)=ExKurtss(SVmnoap);
        Smu(tr)=mean(SVmnoap);
        SSD(tr)=std(SVmnoap);
        
        DSK(tr)=skewns(DVmnoap) ;
        DKur(tr)=ExKurtss(DVmnoap);
        Dmu(tr)=mean(DVmnoap);
        DSD(tr)=std(DVmnoap);
        
        RSK(tr)=skewns(RVmnoap) ;
        RKur(tr)=ExKurtss(RVmnoap);
        Rmu(tr)=mean(RVmnoap);
        RSD(tr)=std(RVmnoap);
        
        RmuMY(tr)=mean(RVmnoapMY);
        DmuMY(tr)=mean(DVmnoapMY);
        SmuMY(tr)=mean(SVmnoapMY);

        
        
        
        SK(tr)=skewns(Vmnoap) ;
        Kur(tr)=ExKurtss(Vmnoap);
        mu(tr)=mean(Vmnoap);
        SD(tr)=std(Vmnoap);
        
        
        figure(100)
        subplot(4,2,1)
        plot(Time,ALLVm)
        hold on
        plot(Time,Vmnoap,'r')
        
        plot(Time,VmnoapMY,'g')
        hold off
        ylim([-90 15])
        
        subplot(4,2,2)
        hist(ALLVm,50)
        xlim([-90 -25])
        
        subplot(4,2,3)
        plot(1:tr,SK(1:tr))
        hold on
        plot(tr,SK(tr),'or')
        ylabel('Skeness')
        
        hold off
        
        subplot(4,2,4)
        plot(1:tr,Kur(1:tr))
        hold on
        plot(tr,Kur(tr),'or')
        ylabel('EKur')
        hold off
        
        
        subplot(4,2,5)
        plot(1:tr,mu(1:tr))
        hold on
        plot(tr,mu(tr),'or')
        
        ylabel('mean')
        hold off
        
        
        
        subplot(4,2,6)
        plot(1:tr ,SD(1:tr))
        hold on
        plot(tr,SD(tr),'or')
        hold off
        ylabel('SD')
        
        
        
        RVtot=[RVtot RVmnoap-mean(RVmnoap)];
        SVtot=[SVtot SVmnoap-mean(SVmnoap)];
        DVtot=[DVtot DVmnoap-mean(DVmnoap)];
        Vtot=[Vtot Vmnoap-mean(Vmnoap)];
        
        VtotNoMeanRed=[VtotNoMeanRed Vmnoap];
        VtotNoMeanRedMY=[VtotNoMeanRedMY VmnoapMY];
        VtotMY=[VtotMY VmnoapMY-mean(VmnoapMY)];
        SVtotMY=[SVtotMY SVmnoapMY-mean(SVmnoapMY)];
        DVtotMY=[DVtotMY DVmnoapMY-mean(DVmnoapMY)];
        RVtotMY=[RVtotMY RVmnoapMY-mean(RVmnoapMY)];
        
        
        %                 pause
    end
    files{ind}.ISI=ISItmp;
    
    % % % % % %     T2=(0:length(DVtot)-1)./sampleRate;
    % % % % % %     tic;    specoutD{ind}=spectralsingleVmNew(T2, DVtot',[0 200]);toc;
    % % % % % %     T2=(0:length(SVtot)-1)./sampleRate;
    % % % % % %     specoutS{ind}=spectralsingleVmNew(T2, SVtot',[0 200]);
    
    files{ind}.TotTreshR=TotTreshR;
    files{ind}.TotTreshS=TotTreshS;
    files{ind}.TotTreshD=TotTreshD;
    files{ind}.TotTreshAll=TotTreshAll;
        
        
    files{ind}.Vtot=Vtot;
    files{ind}.VtotNoMeanRed=VtotNoMeanRed;
    files{ind}.VtotNoMeanRedMY=VtotNoMeanRedMY;
    files{ind}.SVtotMY=SVtotMY;
    files{ind}.DVtotMY=DVtotMY;
    files{ind}.RVtotMY=RVtotMY;
    %                 files{ind}.SVtot=SVtot;
    %                 files{ind}.DVtot=DVtot;
    
    files{ind}.VtotMY=VtotMY;
    
    files{ind}.Vm=Vm;
    files{ind}.Dmu=Dmu;
    files{ind}.DSD=DSD;
    

    files{ind}.Rmu=Rmu;
    files{ind}.RSD=RSD;
    
    TSK(ind)=skewns(Vtot) ;
    TKur(ind)=ExKurtss(Vtot);
    Tmu(ind)=mean(mu);
    TSD(ind)=std(Vtot);
    
    [h,p,kstemp,c] = lillietest(Vmnoap);
    Tkstat(ind)=kstemp; % to make the right comparison
    
    
    
    files{ind}.SVtot=SVtot;
    STSK(ind)=skewns(SVtot) ;
    STKur(ind)=ExKurtss(SVtot);
   
    
    STmu(ind)=mean(Smu);
    STSD(ind)=std(SVtot);
    
    [h,p,kstemp,c] = lillietest(SVmnoap);
    
    STkstat(ind)=kstemp; % to make the right comparison
    
    
    files{ind}.DVtot=DVtot;
    DTSK(ind)=skewns(DVtot) ;
    DTKur(ind)=ExKurtss(DVtot);
    
    DTmu(ind)=mean(Dmu);
    DTSD(ind)=std(DVtot);
    
    
    files{ind}.RVtot=RVtot;
    RTSK(ind)=skewns(RVtot) ;
    RTKur(ind)=ExKurtss(RVtot);   
    RTmu(ind)=mean(Rmu);
    RTSD(ind)=std(RVtot);
    
    DTSDMY(ind)=std(DVtotMY);
    DTmuMY(ind)=mean(DmuMY);
    STmuMY(ind)=mean(SmuMY);
    STSDMY(ind)=std(SVtotMY);
    RTmuMY(ind)=mean(RmuMY);
    RTSDMY(ind)=std(RVtotMY);

    
    [h,p,kstemp,c] = lillietest(DVmnoap);
    
    DTkstat(ind)=kstemp; % to make the right comparison
    
    
    ParDSK(ind)=mean(DSK) ;
    ParDKur(ind)=mean(DKur);
    ParDSD(ind)=mean(DSD);
    ParSSK(ind)=mean(SSK) ;
    ParSKur(ind)=mean(SKur);
    ParSSD(ind)=mean(SSD);
    ParRSK(ind)=mean(RSK) ;
    ParRKur(ind)=mean(RKur);
    ParRSD(ind)=mean(RSD);
    
    nSecD(ind)=length(DVtot)/sampleRate;
    nSecS(ind)=length(SVtot)/sampleRate;
    nSecR(ind)=length(RVtot)/sampleRate;
    
    
    toc;
    
    
    
    FrRALL(ind)=mean(files{ind}.RFr);
    FrSALL(ind)=mean(files{ind}.SFr);
    FrDALL(ind)=mean(files{ind}.DFr);
   
    RCVtmp=files{ind}.RCV;
    RCVtmp(isnan(RCVtmp))=[];
    RCVplot(ind)=mean(RCVtmp);
    
    DCVtmp=files{ind}.DCV;
    DCVtmp(isnan(DCVtmp))=[];
    DCVplot(ind)=mean(DCVtmp);
    
    SCVtmp=files{ind}.SCV;
    SCVtmp(isnan(SCVtmp))=[];
    SCVplot(ind)=mean(SCVtmp);

    figure(1011)
    hold all
    plot(DTmu,DTSD,'ro')
    title('SD Vs Vm')
    
    figure(1012)
    hold all
    plot(DTmu,FrDALL,'ro')
    title('Fr Vs Vm')
    
    figure(1013)
    hold all
    plot(FrDALL,DTSD,'ro')
    title('SD Vs Fr')
    
    
    figure(1011)
    hold all
    % plot(STmu,STSD,'go')
    plot(STmu,STSD,'bo')
    axis tight
    box off
    ylim([0 5])
    
    figure(1012)
    hold all
    plot(STmu,FrSALL,'bo')
    title('Fr Vs Vm')
    
    figure(1013)
    hold all
    % plot(DTmu,DTSD,'bo')
    plot(FrSALL,STSD,'bo')
    title('SD Vs Fr')
    
    
    figure(1011)
    hold all
    plot(DTmu(ind),DTSD(ind),'go')
    
    
    figure(1012)
    hold all
    plot(DTmu(ind),FrDALL(ind),'go')
    title('Fr Vs Vm')
    
    figure(1013)
    hold all
    plot(FrDALL(ind),DTSD(ind),'go')
    title('SD Vs Fr')
    
    
    figure(1011)
    hold all
    % plot(STmu,STSD,'go')
    plot(STmu(ind),STSD(ind),'ko')
    axis tight
    box off
    ylim([0 5])
    
    figure(1012)
    hold all
    plot(STmu(ind),FrSALL(ind),'ko')
    title('Fr Vs Vm')
    
    figure(1013)
    hold all
    % plot(DTmu,DTSD,'bo')
    plot(FrSALL(ind),STSD(ind),'ko')
    title('SD Vs Fr')
    
    
    if isempty(files{ind}.TotTreshAll)
        thS(ind)=0;
    else
        thS(ind)=mean(files{ind}.TotTreshAll);
    end
    
    if isempty(files{ind}.TotTreshAll)
        thD(ind)=0;       
    else
        thD(ind)=mean(files{ind}.TotTreshAll);
    end
    

  
    figure(2011)
    hold all
    plot(DTmuMY,DTSDMY,'ro')
    title('SD Vs Vm')
    
    figure(2012)
    hold all
    plot(DTmuMY,FrDALL,'ro')
    title('Fr Vs Vm')
    
    figure(2013)
    hold all
    plot(FrDALL,DTSDMY,'ro')
    title('SD Vs Fr')
    
    figure(2014)
    hold all
    plot(FrDALL,DCVplot,'ro')
    title('CV Vs Fr')
    
    
    figure(2011)
    hold all
    % plot(STmu,STSD,'go')
    plot(STmuMY,STSDMY,'bo')
    axis tight
    box off
    ylim([0 5])
    
    figure(2012)
    hold all
    plot(STmuMY,FrSALL,'bo')
    title('Fr Vs Vm')
    
    figure(2013)
    hold all
    % plot(DTmu,DTSD,'bo')
    plot(FrSALL,STSDMY,'bo')
    title('SD Vs Fr')
    
     figure(2014)
    hold all
    plot(FrSALL,SCVplot,'bo')
    title('CV Vs Fr')
    
    figure(2011)
    hold all
    plot(DTmuMY(ind),DTSDMY(ind),'go')
    
    
    figure(2012)
    hold all
    plot(DTmuMY(ind),FrDALL(ind),'go')
    title('Fr Vs Vm')
    
    figure(2013)
    hold all
    plot(FrDALL(ind),DTSDMY(ind),'go')
    title('SD Vs Fr')
    
    figure(2014)
    hold all
    plot(FrDALL(ind),DCVplot(ind),'go')
    title('CV Vs Fr')
    
    figure(2011)
    hold all
    % plot(STmu,STSD,'go')
    plot(STmuMY(ind),STSDMY(ind),'ko')
    axis tight
    box off
    ylim([0 5])
    
    figure(2012)
    hold all
    plot(STmuMY(ind),FrSALL(ind),'ko')
    title('Fr Vs Vm')
    
    figure(2013)
    hold all
    % plot(DTmu,DTSD,'bo')
    plot(FrSALL(ind),STSDMY(ind),'ko')
    title('SD Vs Fr')
    
      figure(2014)
    hold all
    plot(FrSALL(ind),SCVplot(ind),'ko')
    title('CV Vs Fr')
    
    
  
%     pause
end

%%
% Then plot:
if LRFlag==2
    indL=1000;
else
    indL=0;
end
mCorrectedDelay=DTmu;

if LRFlag==2
    meanDelayL=DTmu;
    meanEndDelayL2=RTmu; %Response here is going to be END Delay
else
    meanDelayR=DTmu;
    meanEndDelayR2=RTmu; %Response here is going to be END Delay

end
% mCorrectedDelay=DTmu;
% mCorrectedDelay=DTmu;

% for ff=1:89
for ff=1:47

    FrSALL(ff)=mean(files{ff}.SFr);
    FrDALL(ff)=mean(files{ff}.DFr);
end
figure(20+indL)
hold all
% plot(DTmu,DTSK,'bo')
plot(mCorrectedDelay,DTSK,'go')
hold on
x=mCorrectedDelay';
y=(DTSK)';
[Rsq,yCalc]= MyLinearReg(x,y);
plot(x,yCalc,'r')
legend(sprintf('red R^2=%.2f ',Rsq))

title('Skewness Vs Vm')

figure(21+indL)
hold all
% plot(DTmu,DTkstat,'bo')
plot(mCorrectedDelay,DTkstat,'ro')
title('KS Vs Vm')

figure(23+indL)
hold all
% plot(DTmu,DTSD,'bo')
plot(mCorrectedDelay,DTSD,'ro')
title('SD Vs Vm')


x=mCorrectedDelay';
y=(DTSD)';
[Rsq,yCalc]= MyLinearReg(x,y);
plot(x,yCalc,'r')
legend(sprintf('red R^2=%.2f ',Rsq))




figure(112+indL)
hold all
% plot(DTmu,DTSD,'bo')
plot(mCorrectedDelay,FrDALL,'ro')
title('Fr Vs Vm')

figure(113+indL)
hold all
% plot(DTmu,DTSD,'bo')
plot(FrDALL,DTSD,'ro')
title('SD Vs Fr')


figure(32+indL)
hold all
% plot(DTmu,DTSD,'bo')
plot(mCorrectedDelay,DTSD,'ro')
title('SD Vs Vm')


figure(33+indL)
hold all
% plot(DTmu,DTKur,'bo')
plot(mCorrectedDelay,DTKur,'ro')
title('Kurtosis Vs Vm')

figure
plot(mCorrectedDelay(AuditoryFlag(goodCellsNoRepRegular)==0),FrDALL(AuditoryFlag(goodCellsNoRepRegular)==0),'om')
hold on;plot(mCorrectedDelay(AuditoryFlag(goodCellsNoRepRegular)==1),FrDALL(AuditoryFlag(goodCellsNoRepRegular)==1),'xm')

mV=mean(mCorrectedDelay(AuditoryFlag(goodCellsNoRepRegular)==1));SDV=std(mCorrectedDelay(AuditoryFlag(goodCellsNoRepRegular)==1));
mVS=mean(mCorrectedDelay(AuditoryFlag(goodCellsNoRepRegular)==0));SDVS=std(mCorrectedDelay(AuditoryFlag(goodCellsNoRepRegular)==0));
title(sprintf('Aud: mVD=%.1f,sdmVD=%.1f; Tac mVD=%.1f,sdmVD=%.1f',mV,SDV,mVS,SDVS))

figure

plot(mCorrectedDelay(AuditoryFlag(goodCellsNoRepRegular)==0),DTSD(AuditoryFlag(goodCellsNoRepRegular)==0),'om')
hold on;plot(mCorrectedDelay(AuditoryFlag(goodCellsNoRepRegular)==1),DTSD(AuditoryFlag(goodCellsNoRepRegular)==1),'xm')


x=mCorrectedDelay';
y=(DTSD)';
[Rsq,yCalc]= MyLinearReg(x,y);
plot(x,yCalc,'m')
legend(sprintf('red R^2=%.2f ',Rsq))


mV=mean(DTSD(AuditoryFlag(goodCellsNoRepRegular)==1));SDV=std(DTSD(AuditoryFlag(goodCellsNoRepRegular)==1));
mVS=mean(DTSD(AuditoryFlag(goodCellsNoRepRegular)==0));SDVS=std(DTSD(AuditoryFlag(goodCellsNoRepRegular)==0));
title(sprintf('Aud: mSD=%.1f,sdSD=%.1f; Tac mSD=%.1f,sdSD=%.1f',mV,SDV,mVS,SDVS))
box off

% 
% figure
% histogram((mCorrectedDelay(AuditoryFlag(goodCellsNoRepRegular)==1)),-70:4:-30,'facecolor','m','facealpha',.5,'edgecolor','none')
% hold on
% histogram((mCorrectedDelay(AuditoryFlag(goodCellsNoRepRegular)==0)),-70:4:-30,'facecolor','k','facealpha',.5,'edgecolor','none')
% box off

figure
histogram(mCorrectedDelay,-70:4:-30,'facecolor','m','facealpha',.5,'edgecolor','none')
hold on
box off
xlim([-70 -30])

% Plot with respect to threshold and threhsold distribution
difv_plot=mCorrectedDelay;tresh_plot=thD;
mean(thD)
std(thD)
difv_plot(thD==0)=[];tresh_plot(thD==0)=[];
figure
histogram(tresh_plot,-40:0.5:-30,'facecolor','m','facealpha',.5,'edgecolor','none')

figure
histogram(tresh_plot-difv_plot,6:2:30,'facecolor','m','facealpha',.5,'edgecolor','none')






%% Figure S1C: CV Vs Rates 

load('FrCVWholeDelay.mat')
figure(300)
plot(EtotFrD,EtotCVD,'ob')
hold on
plot(FrDALL,CVD,'xg')
plot(EtotFrS,EtotCVS,'or')
hold on
plot(FrSALL,CVS,'xc')
box off
legend('Delay,Probes','Delay,WholeCell','Sample,Probes','Sample,WholeCell')
xlim([1 60])
ylim([0 3])
