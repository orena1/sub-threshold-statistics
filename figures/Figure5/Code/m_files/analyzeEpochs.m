function analyzeEpochs(iwdata,params, name, ctrial)

% derived from analyzeVmfluct5

if nargin<4
    ctrial=1;
    if nargin<3
        name=[];
    end;
end;


if isempty(params)
    
    params.touchPSPdur=0.1; % in seconds
    params.prtl=5; % 5th percentile.
    params.ampth=2.5;
    params.md=0.025; % when searching for peaks, minimal separation distance
    params.plot_ratio=1;
    params.vrange=[-70 -20]; % plot vm range
    params.usemask=0;
    params.range=[]; % range is the beginning and ending of counting window. if empty you count the whole trace
    params.laser=[];
end;
% 
% if ~isfield (params, 'Vm')
%     params.Vm=1; % when params.Vm=1, it is membrane potential analysis instead of spike analysis.
% end;

% if ~isfield(params, 'plotexamples')
%     params.plotexamples=5;
% end;

touchPSPdur=params.touchPSPdur;
prtl=params.prtl;
ampth=params.ampth;
md=params.md;
plot_ratio=params.plot_ratio;
vrange=params.vrange;
usemask=params.usemask;

% this program will analyze a more genuine form of data structure.
% iwdata2 =
%
%         cellname: 'JY0520AAAC'
%          mainwid: 1
%           allwid: [0 1]
%        trialnums: [1x75 double]
%                k: 75
%                t: [1x5000 double]
%              tvm: [1x52000 double]
%              Vth: [1x1 struct]
%            S_ctk: [4-D double]
%     featuresName: {1x11 cell}
%               Vm: [2x5000x75 double]
%            Vmorg: [52000x75 double]
%              Spk: [5000x75 double]
%  wdata.varNames={ 'thetaAtBase'  'amplitude'  'phase'    'setpoint'  'deltaKappa'    'moment'    'FaxialAdj', 'touch_onset', 'touch_offset', 'licks'};
%  iwdata2.featuresName
% will also store the selections.
% r is the mask

% modified 8.27.2014, use amplitude alone to find out "nonwhisking"
% periods.

% prtl: percentile selection, e.g., 10th, this is the chosen baseline
% wid: whisker id
% prefph: preferred phase. if over pi, it means going through all whisking
% phases

% touchPSPdur about 200 ms, where to search the maximum depolarization

% ampth: amplitude threshold, where whisking is considered to be whisking

% md: min distance separation between peaks. in ms

% plot_ratio, if 0.25, means 25% of all figures will be plotted.

% analyze the amplitude of Vm fluctuations between truly-spontaneous,
% whisking, and touch periods. The baseline is the resting potential.

% nontrials are those that are not included, these could be stimulation
% trials, unstable trials, etc.
% every three trials, there will be an estimation of the baseline

attendedcells={'JY1036'  % contains a large number of spiklets
    };

Fs=round(1/median(diff(iwdata.tvm))); % sampling rate of Vm

% w_params=whiskdecomposej(wpos, tpos, twhiskorg);

trialnums=iwdata.trialnums;

Vm=iwdata.Vmorg;
tvm=iwdata.tvm; tvm=tvm';
Spk=iwdata.Spkorg;

if ~isempty(params.laser)
    switch params.laser
        
        case 'b'
            laser=iwdata.opto(:, :, 1);
        case 'o'
            if ndims(iwdata.opto)>2
                laser=iwdata.opto(:, :, 2);
            else
                laser=iwdata.opto;
            end;
        otherwise
            laser=zeros(size(Vm));
    end;
else
    
    laser=zeros(size(Vm));
end;

indlimit=find(tvm<=5);

Vm=Vm(indlimit, :);
Spk=Spk(indlimit, :);
laser=laser(indlimit, :);
tvm=tvm(indlimit);
Spk=iwdata.Spkorg;

baseline=zeros(1, size(Vm, 2));

PSP_nw_all=[];
PSP_wh_all=[];

% PSP_TOUCH=cell(1, length(trialnums));
% PSP_touchtime=cell(1, length(trialnums));
% PSP_touch_all=cell(1, length(iwdata.allwid));

% new parameter: dtt means distance-to-threshold
% dtt_nw_all=[];
% dtt_wh_all=[];
% dtt_TOUCH_all=cell(1, length(iwdata.allwid));

Vmtransit=[];% collect nonwhisking-to-whisking transition, 0.2s pre and 0.5s post, if touch occurs, filled the data with nan
Vmtransit2=[];
Spktransit=[];
Whisktransit=[];
Whiskpostransit=[];

transit_index=[]; % trial number of transit, time of the transit. 

% now go through each trial

Vmph={}; % first row is alway whisking, second row is always Vm (raw, with spikes)
Spkph={};
phVm={};

% collect all data.
Epochs.twhiskVm=[];
Epochs.whiskVm=[];
Epochs.whiskSpk=[];
Epochs.whiskamp=[];
Epochs.whiskangle=[];
Epochs.whiskphase=[];

Epochs.whisk_index=[]; % Here, trial number and onset of the whisk epoch should be written down
Epochs.nonwhisk_index=[]; % same, trial number and onset of non-whisk epoch. 

Epochs.tnonwhiskVm=[];
Epochs.nonwhiskVm=[];
Epochs.nonwhiskSpk=[];
Epochs.nonwhiskamp=[];

h2=waitbar(0, ['Analyzing...' iwdata.cellname]);
ctrial=1;
for i=ctrial:length(trialnums)
    
    waitbar(i/length(trialnums), h2)
    itrialnum=trialnums(i);
    ilaser=laser(:, i);
    iVmorg=Vm(:, i);
    iVm=sgolayfilt(iVmorg, 3, 21);
    iSpk=Spk(:, i);
    iSpktime=iwdata.tvm(find(iwdata.Spkorg(:, i))); % spike time
%     iVth=iwdata.Vth.threshold{i}; % threshold at spike time
%     if ~isempty(iVth)
%         spkthfixed=prctile(iVth, 25); % a fixed spike threshold value.
%     elseif ~isempty(cell2mat(iwdata.Vth.threshold(max(i-5, 1):min(i+5, length(trialnums)))))&&length(cell2mat(iwdata.Vth.threshold(max(i-5, 1):min(i+5, length(trialnums)))))>1
%         spkthfixed=prctile(cell2mat(iwdata.Vth.threshold(max(i-5, 1):min(i+5, length(trialnums)))), 25);
%     elseif ~isempty(cell2mat(iwdata.Vth.threshold(max(i-10, 1):min(i+3, length(trialnums)))))&& length(cell2mat(iwdata.Vth.threshold(max(i-10, 1):min(i+10, length(trialnums)))))>1
%         spkthfixed=prctile(cell2mat(iwdata.Vth.threshold(max(i-10, 1):min(i+10, length(trialnums)))), 25);
%     elseif ~isempty(cell2mat(iwdata.Vth.threshold)) && length(cell2mat(iwdata.Vth.threshold))>1
%         spkthfixed=prctile(cell2mat(iwdata.Vth.threshold), 25);
%     else
%         spkthfixed=-40; % if no spike th can be found, use -44 mV as a fixed value.
%     end;
%     
    %     iVmorg=sgolayfilt(iVmorg, 3, 2
    if ~isempty(iwdata.Vthparams)
        
        if isempty(iSpktime)
            iVmnoap=sgolayfilt(iVmorg, 3, 21);
        else
            iVmnoap=sgolayfilt(removeAPnewP(iVmorg,Fs,0.33,  iwdata.Vthparams, 10000), 3, 21); % remove AP
        end;
        
    else
        iVmnoap=iVmorg;
    end;
    
    if ~isempty(find(ismember(iwdata.cellname(1:6), attendedcells)))
        iVmnoap=medfilt1(iVmnoap, 51);
    end;
    
    Vm_clust=sgolayfilt(iVmnoap, 3, 51);
    baseline(i)=prctile(Vm_clust, prtl);
    
    if length(size(iwdata.S_ctk))==4
        %           for two whiskers:
        %           touchoff              5000x2
        %           touchon               5000x2
        % iwdata.t(find(touchon(:, wid)))--> touch time on whisker wid
        touchon=squeeze(iwdata.S_ctk(8, :, i, :)); %
        touchoff=squeeze(iwdata.S_ctk(9, :, i, :));
    else
        touchon=squeeze(iwdata.S_ctk(8, :, i)); % is a single vector
        touchoff=squeeze(iwdata.S_ctk(9, :, i));
        
        if size(touchon, 1)<size(touchon, 2)
            touchon=touchon';
            touchoff=touchoff';
        end;
    end
    
    % licking time, need them so that we won't include them for whiksing or
    % nonwhisking
    licktime=squeeze(iwdata.S_ctk(10, :, :, 1));
    %% touch PSPs delta Vm: see how much depolarization touchPSPdur sec after touch, from baseline
    
%     PSP_touchtime{i}=cell(1, length(iwdata.allwid));
%     PSP_TOUCH{i}=cell(1, length(iwdata.allwid));
%     dtt_TOUCH{i}=cell(1, length(iwdata.allwid));
%     %     collection of Vm and Spk
    %     touchVm=[];
    %     touchSpk=[];
%     
    for iw=1:size(touchon, 2) % if only one whisker, size(touchon, 2)=1
        touch_mask_compute{iw}=sparse(size(tvm, 1), size(tvm, 2));% this mask will be used to compute touch PSP
        if ~isempty(find(touchon(:, iw))) % there is touch on this whisker
            % construct a touch mask based on whisker iw
            touchtime{iw}=iwdata.t(find(touchon(:, iw))); % time of touch, in sec
            touchofftime{iw}=iwdata.t(find(touchoff(:, iw)));
            for k=1:length(touchtime{iw})
                touch_mask_compute{iw}(tvm>touchtime{iw}(k) & tvm<touchtime{iw}(k)+touchPSPdur)=1;
            end;
            touch_eps=dissectmat(touch_mask_compute{iw}, 10);
            %
            %                 for ite=1:size(touch_eps)
            %
            %                     tPSP=tvm(touch_eps(ite, 1):touch_eps(ite, 2));
            %                     vmPSP=iVmnoap(touch_eps(ite, 1):touch_eps(ite, 2));
            %                     vmPSPorg=iVm(touch_eps(ite, 1):touch_eps(ite, 2));
            %                     spkPSP=iSpk(touch_eps(ite, 1):touch_eps(ite, 2));
            %
            %                     if iwdata.allwid(iw)==iwdata.mainwid
            %                         % collect these for storage.
            %                         Epochs.touchVm=[Epochs.touchVm {vmPSP}];
            %                         Epochs.touchSpk=[Epochs.touchSpk {spkPSP}];
            %                     end;
            %                     % if spike occurs, PSP peak is the spike threshold
            %                     % otherwise, PSP peak is the real peak.
            %
            %                     if any((tPSP(1)-iSpktime).*(tPSP(end)-iSpktime)<0) % spike osccurs within tPSP
            %                         % here are the first spike within tPSP
            %                         spktime_touch=iSpktime(iSpktime>=tPSP(1)&iSpktime<=tPSP(end));
            %
            %                         Vth_first=iVth(find(iSpktime>tPSP(1), 1, 'first'));
            %                         PSPpeakmax=min(iVth(iSpktime>=tPSP(1) & iSpktime<=tPSP(end))); % first spikes
            %                         [PSP_touch_peaks, touch_peaks_locs]=max(vmPSP); % only choose one max depolarization per touch.
            %
            %                         tPSPpeak=tPSP(touch_peaks_locs);
            %                         dtt_TOUCH_ite=min(0, -Vth_first+PSP_touch_peaks);% any point exceeding the threshold is 0
            %
            %                         % also update the peak time.
            %                         PSPspkind=find(PSP_touch_peaks>=PSPpeakmax);
            %                         if ~isempty(PSPspkind)
            %                             PSP_touch_peaks(PSP_touch_peaks>=PSPpeakmax)=PSPpeakmax;
            %                             for ip=1:length(PSPspkind)
            %                                 if any(find(spktime_touch<tPSPpeak(PSPspkind(ip)), 1, 'first'))
            %                                     tPSPpeak(PSPspkind(ip))=spktime_touch(find(spktime_touch<tPSPpeak(PSPspkind(ip)), 1, 'first'));
            %                                 end;
            %                             end;
            %                         end;
            
            %                     figure (50); clf
            %                     plot(tPSP, vmPSPorg, 'k');
            %                     hold on
            %                     plot(tPSP, vmPSP, 'color', [0.5 0.5 0.5]);
            %                     line([tPSPpeak, tPSPpeak], [PSP_touch_peaks PSP_touch_peaks+dtt_TOUCH_ite], 'color', 'r', 'linewidth', 2)
            %
            
            % all touches on this trial
            %                         PSP_touchtime{i}{iw}=[PSP_touchtime{i}{iw}; tPSPpeak'];
            %                         PSP_TOUCH{i}{iw}=[PSP_TOUCH{i}{iw}; PSP_touch_peaks-baseline(i)];
            %                         PSP_touch_all{iw}=[PSP_touch_all{iw}; PSP_touch_peaks-baseline(i)];
            %
            %                         dtt_TOUCH{i}{iw}=[dtt_TOUCH{i}{iw}; dtt_TOUCH_ite];
            %                         dtt_TOUCH_all{iw}=[dtt_TOUCH_all{iw};  dtt_TOUCH_ite];
            %
            %                     else
            %                         [PSP_touch_peaks, touch_peaks_locs]=max(vmPSP);
            %                         tPSPpeak=tPSP(touch_peaks_locs);
            %
            %                         dtt_TOUCH_ite=-spkthfixed+PSP_touch_peaks;
            %
            %                         PSP_touchtime{i}{iw}=[PSP_touchtime{i}{iw}; tPSPpeak'];
            %                         PSP_TOUCH{i}{iw}=[PSP_TOUCH{i}{iw}; PSP_touch_peaks-baseline(i)];
            %                         PSP_touch_all{iw}=[PSP_touch_all{iw}; PSP_touch_peaks-baseline(i)];
            %
            %                         dtt_TOUCH{i}{iw}=[dtt_TOUCH{i}{iw}; dtt_TOUCH_ite];
            %                         dtt_TOUCH_all{iw}=[dtt_TOUCH_all{iw};  dtt_TOUCH_ite];
            %
            %                     end;
            %
            %            end;
        end;
        %                         figure(66); clf
        %                         plot(tvm, iVmnoap);hold on
        %                         plot(tvm(find(touch_mask_compute{iw})), iVmnoap(find(touch_mask_compute{iw})), 'k.');
        %                         line([min(tvm) max(tvm)], [baseline(i) baseline(i)], 'color', 'k')
        %                         if ~isempty(PSP_touchtime{i}{iw})
        %                         plot(PSP_touchtime{i}{iw}, PSP_TOUCH{i}{iw}+baseline(i),'ro', 'markersize', 8, 'linewidth', 2.5)
        %                         end;
        %                         pause
        %         touchofftime{iw}=iwdata.t(find(touchoff(:, iw))); %
    end;
%     
    
    %% making masks from touch and get whisking information
    % masks for touch
    % if any nonwhisking or whisking episodes happens in touch_mask, it will not be counted.
    
    touch_mask=[];
    if ~isfield(params, 'touchmask')
        params.touchmask=1;
    end;
    
    if params.touchmask
        for iw=1:size(touchon, 2)
            touch_mask=[touch_mask; -0.01+iwdata.t(find(touchon(:, iw)))' 0.1+iwdata.t(find(touchoff(:, iw)))']; % in seconds
        end;
    end;
    
    % next, transform touch_mask into index that matches with tvm
    touch_mask_ind=sparse(size(tvm, 1), size(tvm, 2));
    if ~isempty(touch_mask)
        for it=1:size(touch_mask, 1)
            touch_mask_ind(tvm>=touch_mask(it, 1) & tvm<=touch_mask(it, 2))=1; % all touch periods are 1, else are 0
        end;
    end;
    
    twhisk=iwdata.t;
    licks=[];
    if length(iwdata.allwid)>1
        whisktheta=         squeeze(iwdata.S_ctk(1, :, i, iwdata.mainwid+1));           % whisking angle
        whiskthetafilt=     squeeze(iwdata.S_ctk(5, :, i, iwdata.mainwid+1));
        whiskamp=           squeeze(iwdata.S_ctk(2, :, i, iwdata.mainwid+1));           % whisking amplitude
        whiskphase=         squeeze(iwdata.S_ctk(3, :, i, iwdata.mainwid+1));           % whisking phase
        licks=              iwdata.t(find(squeeze(iwdata.S_ctk(10, :, i, 1))));               % licks in seconds
    else
        whisktheta=         squeeze(iwdata.S_ctk(1, :, i));                                      % whisking angle
        whiskthetafilt=     squeeze(iwdata.S_ctk(5, :, i));
        whiskamp=           squeeze(iwdata.S_ctk(2, :, i));                                      % whisking amplitude
        whiskphase=         squeeze(iwdata.S_ctk(3, :, i));                                      % whisking phase
        licks=              iwdata.t(find(squeeze(iwdata.S_ctk(10, :, i))));                  % licks in seconds
    end;
    
    % licks masks are here
    lick_mask_ind=sparse(size(tvm, 1), size(tvm, 2));
    
    if ~isfield(params, 'lickmask')
        params.lickmask=1;
    end;
    
    if params.lickmask==1
        if ~isempty(licks)
            for it=1:numel(licks)
                lick_mask_ind(tvm>=licks(it)-0.1 & tvm<=licks(it)+0.10)=1; % all touch periods are 1, else are 0
            end;
        end;
    end;
    
    % another window for restriction the analysis window.
    count_mask_ind=sparse(size(tvm, 1), size(tvm, 2));
        
    if ~isempty(params.laser)
        
        laseron=find(ilaser>1);
        
        onsets=laseron([1; 1+find(diff(laseron)>1)]);
        laserstimfreq=10000/median(diff(onsets));
        
        if laserstimfreq>=20
            
            count_mask_ind([1:laseron(1)+params.delay*10000 laseron(end):end])=1;
            
        else
            count_mask_ind(find(ilaser<1))=1;
        end;
    end;
    
    if ~isempty(params.range)
        count_mask_ind(tvm<params.range(1) | tvm>params.range(2))=1;
    end;
    
    %% nonwhisking periods
    % Fill this one out 
  
    nw_mask_ind_new=sparse(size(tvm, 1), size(tvm, 2));
    %     if ~usemask
    % passing amplitude
    nw_periods = denovodetectWhiskingEpochs_general(iwdata, i, 1.25, 0.1, 0); % remember to change i+1 back to i!!!
    % touch should not be within "non whisking periods" or 150 ms near it
    nw_mask_ind=sparse(size(tvm, 1), size(tvm, 2));
    if ~isempty(nw_periods{1})
        % generate non-whisking mask
        for inw=1:size(nw_periods{1}, 1)
            nw_mask_ind(tvm>=nw_periods{1}(inw, 1) & tvm<=nw_periods{1}(inw, 2))=1; % all touch periods are 1, else are 0
        end
    end;
    
    % get rid of the lick and touch overlaps
    nw_mask_ind(find(touch_mask_ind))=0;
    nw_mask_ind(find(lick_mask_ind))=0;
    nw_mask_ind(find(count_mask_ind))=0;
    
    % examine to get rid of short periods.
    if ~isempty(find(nw_mask_ind))
        nwinds=find(nw_mask_ind);
        diff_nw=diff(nwinds);
        
        nw_start=nwinds([1 1+find(diff_nw>1)]);
        nw_end=[nwinds(find(diff_nw>1)) nwinds(end)];
        nw_dur=nw_end-nw_start;
        ind_realnw=find(nw_dur>0.1*Fs);
        nw_start=nw_start(ind_realnw);
        nw_end=nw_end(ind_realnw);
    else
        nw_start=[];
        nw_end=[];
    end;
   
    % Epochs.nonwhisk_index=[]; % same, trial number and onset of non-whisk epoch. 
    if ~isempty(nw_start)
        for kk=1:numel(nw_start)
            nw_mask_ind_new(tvm>=tvm(nw_start(kk)) & tvm<=tvm(nw_end(kk)))=1;
        end;
    end;
    
    if ~isempty(nw_mask_ind_new)
        nw_mask_ind_new=mergegaps(nw_mask_ind_new, .1*Fs, 2); % merge gaps less than 100 ms
    end;
    
    
    %% Here, extract transition between nonwhisking and whisking,
    % requirement is within 250 ms after the end of nonwhisking, there should
    % be no touch
    % if there is touch after 250 ms, make the data afterwards NaN
    
    % Vmtransit=[];% collect nonwhisking-to-whisking transition, 0.2s pre and 0.5s post, if touch occurs, filled the data with nan
    % Vmtransit2=[];
    % Spktransit=[];
    % Whisktransit=[];
    % Whiskpostransit=[];
    
    % for current trial.
    PSP_NW{i}=[];
    PSPNW_time{i}=[];
    
    dtt_NW{i}=[];
    
    for inw=1:length(nw_start)
        
        if tvm(nw_end(inw))-tvm(nw_start(inw))> 0.2 && tvm(nw_end(inw))>=0.201 && tvm(nw_end(inw))+0.501<=max(tvm) % non-whisking longer than 200 ms, begining and end within the data length
            
            
            ind_vm=nw_end(inw);
            
            if isempty(find(count_mask_ind( ind_vm-0.2*Fs : ind_vm+0.5*Fs)==1))
            
            [~, ind_whisk]=min(abs(iwdata.t-tvm(nw_end(inw))));
            
            Vmprime=iVmorg;             Vmprime([find(touch_mask_ind) find(lick_mask_ind)])=NaN;
            Vmnoapprime=iVmnoap;        Vmnoapprime([find(touch_mask_ind) find(lick_mask_ind)])=NaN;
            Spkprime=iSpk;              Spkprime([find(touch_mask_ind) find(lick_mask_ind)])=NaN;
            
            vm_transit_inw=             Vmprime(ind_vm-0.2*Fs : ind_vm+0.5*Fs);
            vm_transit_inw2=            Vmnoapprime(ind_vm-0.2*Fs : ind_vm+0.5*Fs);
            spk_transit_inw=            Spkprime(ind_vm-0.2*Fs : ind_vm+0.5*Fs);
            tvm_transit_inw=            tvm(ind_vm-0.2*Fs : ind_vm+0.5*Fs);
            
            whisk_transit_inw=          whiskamp(ind_whisk-0.2*1000 : ind_whisk+0.5*1000);
            whiskpos_transit_inw=       whiskthetafilt(ind_whisk-0.2*1000 : ind_whisk+0.5*1000);
            twhisk_transit_inw=         iwdata.t(ind_whisk-0.2*1000 : ind_whisk+0.5*1000);
            
            Vmtransit=[Vmtransit vm_transit_inw];% collect nonwhisking-to-whisking transition, 0.2s pre and 0.5s post, if touch occurs, filled the data with nan
            Vmtransit2=[Vmtransit2 vm_transit_inw2];
            Spktransit=[Spktransit spk_transit_inw];
            Whisktransit=[Whisktransit whisk_transit_inw'];
            Whiskpostransit=[Whiskpostransit whiskpos_transit_inw'];
            transit_index=[transit_index;  itrialnum  tvm(nw_end(inw))]; % trial number of transit, time of the transit. 
            end;
        end;
        
        tvm_inw=        tvm(nw_start(inw):nw_end(inw));
        vm_orginw=      iVmorg(nw_start(inw):nw_end(inw));
        vm_inw=         iVmnoap(nw_start(inw):nw_end(inw));
        % now based on
        %%
        % **nw_mask_ind_new, nw_start, nw_end** , find out PSP
        % distribution
        Epochs.tnonwhiskVm=     [Epochs.tnonwhiskVm {tvm_inw}];
        Epochs.nonwhiskVm=     [Epochs.nonwhiskVm     {vm_orginw}];
        Epochs.nonwhiskSpk=    [Epochs.nonwhiskSpk    {iSpk(nw_start(inw):nw_end(inw))}];
        Epochs.nonwhiskamp=    [Epochs.nonwhiskamp    {whiskamp(iwdata.t>=tvm(nw_start(inw))&iwdata.t<=tvm(nw_end(inw)))}];
        Epochs.nonwhisk_index=[Epochs.nonwhisk_index; itrialnum tvm(nw_start(inw))]; % given an index. 
        
        %         if params.Vm
        %
        %             if any((tvm_inw(1)-iSpktime).*(tvm_inw(end)-iSpktime)<0)
%                 
%                 spktime_inw=iSpktime(iSpktime>=tvm_inw(1)&iSpktime<=tvm_inw(end));
%                 spkth_inw=iVth(iSpktime>=tvm_inw(1)&iSpktime<=tvm_inw(end));
%                 
%                 [PSP_nw_peaks, nw_peaks_locs]=findpeaks(vm_inw, 'minpeakdistance', md);
%                 
%                 dtt_nw_inw=-min(spkth_inw)+PSP_nw_peaks; % distance to spike threshold is here32
%                 dtt_nw_inw(dtt_nw_inw>0)=0;
%                 
%                 if ~isempty(PSP_nw_peaks>min(spkth_inw))
%                     PSP_nw_peaks(PSP_nw_peaks>min(spkth_inw))=min(spkth_inw);
%                 end;
%                 
%                 tPSPpeak=tvm_inw(nw_peaks_locs);
%                 
%                 % also update the peak time.
%                 PSPspkind=find(PSP_nw_peaks>min(spkth_inw));
%                 if ~isempty(PSPspkind)
%                     PSP_nw_peaks(PSP_nw_peaks>min(spkth_inw))=min(spkth_inw);
%                     for ip=1:length(PSPspkind)
%                         if any(find(spktime_inw<tPSPpeak(PSPspkind(ip)), 1, 'first'))
%                             tPSPpeak(PSPspkind(ip))=spktime_inw(find(spktime_inw<tPSPpeak(PSPspkind(ip)), 1, 'first'));
%                         end;
%                     end;
%                 end;
%                 
%                 PSP_NW{i}=[PSP_NW{i} ;PSP_nw_peaks-baseline(i)];
%                 dtt_NW{i}=[dtt_NW{i}; dtt_nw_inw];
%                 PSPNW_time{i}=[PSPNW_time{i} tPSPpeak]; % peak time in this trial, for plotting purpose
%                 PSP_nw_all=[PSP_nw_all; PSP_nw_peaks-baseline(i)];
%                 dtt_nw_all=[dtt_nw_all; dtt_nw_inw];
%             else
%                 
%                 [PSP_nw_peaks, nw_peaks_locs]=findpeaks(vm_inw, 'minpeakdistance', md);
%                 dtt_nw_inw=-spkthfixed+PSP_nw_peaks; dtt_nw_inw(dtt_nw_inw>0)=0;
%                 PSP_NW{i}=[PSP_NW{i}; PSP_nw_peaks-baseline(i)];
%                 
%                 PSPNW_time{i}=[PSPNW_time{i} tvm_inw(nw_peaks_locs)]; % peak time in this trial, for plotting purpose
%                 PSP_nw_all=[PSP_nw_all; PSP_nw_peaks-baseline(i)];
%                 dtt_nw_all=[dtt_nw_all; dtt_nw_inw];
%             end;
%         end;
    end;
    %---------------------------------------
    %     figure;
    %     subplot(2, 1, 1)
    %     plot(iwdata.t, squeeze(iwdata.S_ctk(2, :, i, 2)))
    %     hold on
    %     if ~isempty(find(nw_mask_ind_new, 1))
    %         plot(tvm(find(nw_mask_ind_new)), 10, 'r.')
    %         if params.usemask
    %             plot(iwdata.t(find(allnw)), 10, 'mo')
    %         end;
    %     end;
    %     line([0 5], [1.25 1.25], 'color', 'r')
    %     axis tight
    %
    %
    %     subplot(2, 1, 2)
    %     plot(tvm, iVmorg, 'k');
    %     hold on
    %     plot(PSP_time{i}, PSP_NW{i}+baseline(i), 'ro', 'linewidth', 2);
    %     line([min(tvm) max(tvm)], [baseline(i) baseline(i)], 'color', 'b')
    %     axis tight
    %----------------------------------------
    
    %% whisking periods, extraction
    
   % if ~usemask % won't use masks. de novo
        indnearph0=find(abs(whiskphase)<0.1);
        indnearph0=indnearph0([1 1+find(diff(indnearph0)>1)]); % this is supposely the peak of protraction
        % remove small or pseudo whisking peaks
        indnearph0(whiskamp(indnearph0)<ampth)=[];
        indnearph0(iwdata.t(indnearph0)<0.05 & iwdata.t(indnearph0)>4.95)=[];
        whisk_mask_ind=sparse(size(tvm, 1), size(tvm, 2));
        
        % build whisking masks
        diffph=[0 diff(whiskphase)];
        cyc_beg_tvm=[];
        cyc_end_tvm=[];
        cycle_minphase=find(diffph<-4); % this is the beginning af all the whisking cycles
        
        cycle_maxphase=cycle_minphase-1;
        cycle_minphase=cycle_minphase(1:end-1);
        cycle_maxphase=cycle_maxphase(2:end);
        
        for iph0=1:length(indnearph0)
            % find out which cycl_minphase proceeds the ph0 point
            ind_beg=find(whiskphase>whiskphase(indnearph0(iph0))&iwdata.t<iwdata.t(indnearph0(iph0))&whiskphase>2, 1, 'last')+1;
            ind_end=find(whiskphase<whiskphase(indnearph0(iph0))&iwdata.t>iwdata.t(indnearph0(iph0))&whiskphase<-2, 1, 'first')-1;
            if ~isempty(ind_beg) && ~isempty(ind_end)
                if ind_end-ind_beg>25
                    whisk_mask_ind(tvm>=iwdata.t(ind_beg)& tvm<=iwdata.t(ind_end))=1;
                end;
            end;
        end;
        if ~isempty(whisk_mask_ind)
            whisk_mask_ind=mergegaps(whisk_mask_ind, .1*Fs, 2); % gaps less than 50 ms will be filled.
        end;
%                 figure;
%                 plot(iwdata.t, whiskphase, 'k.')
%         
%                 hold on
%                 plot(iwdata.t, detrend(whisktheta)/std(whisktheta), 'r')
%                 if ~isempty(find(whisk_mask_ind, 1))
%                     plot(tvm(find(whisk_mask_ind)), 0, 'b.')
%                 end
        
        whisk_mask_ind(find(touch_mask_ind))=0;
        whisk_mask_ind(find(lick_mask_ind))=0;
        whisk_mask_ind(find(count_mask_ind))=0;
        
        if ~isempty(find(whisk_mask_ind))
            [~, whisk_mask_ind]=dissectmat(whisk_mask_ind, .1*Fs); % cycle or bout length has to be larger than 50 ms
        end
   
    
    if ~isempty(find(whisk_mask_ind))
        whisk_mask_ind=mergegaps(whisk_mask_ind, .1*Fs, 2); % merge gaps less han 100 ms.
        [wh_episodes, whisk_mask_ind]=dissectmat(whisk_mask_ind, .1*Fs);
        
        for ie=1:size(wh_episodes)
            
            tvm_iwh=           tvm(wh_episodes(ie, 1):wh_episodes(ie, 2));
            vmorg_iwh=         iVmorg(wh_episodes(ie, 1):wh_episodes(ie, 2));
            vm_iwh=            iVmnoap(wh_episodes(ie, 1):wh_episodes(ie, 2));
            
            % The following was set up at the beginning
            % whiskVm=[];
            % whiskSpk=[];{iSpk(nw_start(inw):nw_end(inw))}
            % whiskamp=[];
            % pile this up for future analysis
            Epochs.twhiskVm=        [Epochs.twhiskVm {tvm_iwh}];
            Epochs.whiskVm=        [Epochs.whiskVm      {vmorg_iwh}];
            Epochs.whiskSpk=       [Epochs.whiskSpk     {iSpk(wh_episodes(ie, 1):wh_episodes(ie, 2))}];
            Epochs.whiskamp=       [Epochs.whiskamp     {whiskamp(iwdata.t>=tvm(wh_episodes(ie, 1))&iwdata.t<=tvm(wh_episodes(ie, 2)))}];
            Epochs.whiskangle=     [Epochs.whiskangle   {whisktheta(iwdata.t>=tvm(wh_episodes(ie, 1))&iwdata.t<=tvm(wh_episodes(ie, 2)))}];
            Epochs.whiskphase=     [Epochs.whiskphase   {whiskphase(iwdata.t>=tvm(wh_episodes(ie, 1))&iwdata.t<=tvm(wh_episodes(ie, 2)))}];
            Epochs.whisk_index=     [Epochs.whisk_index; itrialnum  tvm_iwh(1)]; % given an index.
            
            
            % testing purpose
            %             figure(33); clf (33); plot(tvm_iwh, vmorg_iwh)
            %
%             if params.Vm
%                 if any((tvm_iwh(1)-iSpktime).*(tvm_iwh(end)-iSpktime)<0)
%                     spktime_iwh=    iSpktime(iSpktime>=tvm_iwh(1)&iSpktime<=tvm_iwh(end));
%                     spkth_iwh=      iVth(iSpktime>=tvm_iwh(1)&iSpktime<=tvm_iwh(end));
%                     
%                     if length(vm_iwh)>md
%                         [PSP_wh_peaks, wh_peaks_locs]=findpeaks(vm_iwh, 'minpeakdistance', md);
%                     else
%                         [PSP_wh_peaks, wh_peaks_locs]=max(vm_iwh);
%                     end;
%                     
%                     dtt_wh_iwh=PSP_wh_peaks-min(spkth_iwh); % distance to spike threshold is here32
%                     dtt_wh_iwh(dtt_wh_iwh>0)=0;
%                     
%                     tPSPpeak=tvm_iwh(wh_peaks_locs);
%                     
%                     % also update the peak time.
%                     PSPspkind=find(PSP_wh_peaks>min(spkth_iwh));
%                     if ~isempty(PSPspkind)
%                         PSP_wh_peaks(PSP_wh_peaks>min(spkth_iwh))=min(spkth_iwh);
%                         for ip=1:length(PSPspkind)
%                             if any(find(spktime_iwh<tPSPpeak(PSPspkind(ip)), 1, 'first'))
%                                 tPSPpeak(PSPspkind(ip))=spktime_iwh(find(spktime_iwh<tPSPpeak(PSPspkind(ip)), 1, 'first'));
%                             end
%                         end;
%                     end;
%                     
%                     PSP_WH{i}=[PSP_WH{i} ;PSP_wh_peaks-baseline(i)];
%                     PSPwh_time{i}=[PSPwh_time{i}  tPSPpeak]; % peak time in this trial, for plotting purpose
%                     PSP_wh_all=[PSP_wh_all; PSP_wh_peaks-baseline(i)];
%                     dtt_WH{i}=[dtt_WH{i}; dtt_wh_iwh];
%                     dtt_wh_all=[dtt_wh_all; dtt_wh_iwh];
%                     
%                 else
%                     if length(vm_iwh)>md
%                         [PSP_wh_peaks, wh_peaks_locs]=findpeaks(vm_iwh, 'minpeakdistance', md);
%                     else
%                         [PSP_wh_peaks, wh_peaks_locs]=max(vm_iwh);
%                     end;
%                     
%                     dtt_wh_iwh=-spkthfixed+PSP_wh_peaks; dtt_wh_iwh(dtt_wh_iwh>0)=0;
%                     PSP_WH{i}=[PSP_WH{i} ;PSP_wh_peaks-baseline(i)];
%                     PSPwh_time{i}=[PSPwh_time{i} tvm_iwh(wh_peaks_locs)]; % peak time in this trial, for plotting purpose
%                     PSP_wh_all=[PSP_wh_all; PSP_wh_peaks-baseline(i)];
%                     dtt_WH{i}=[dtt_WH{i}; dtt_wh_iwh];
%                     dtt_wh_all=[dtt_wh_all; dtt_wh_iwh];
%                 end;
%             end;
        end;
    end;
    
    % ---------------------------------------
    % plot the result for saving purpose.
%     
%     if params.Vm==1  ||  i<params.plotexamples
 
if rand>(1-plot_ratio)
    
    hf(i)=figure;
    set(hf(i), 'units', 'centimeter', 'position', [1 2 12 15],'paperpositionmode', 'auto', 'color', 'w')
    
    ha1=subplot(2, 1, 1);
    set(gca, 'nextplot', 'add', 'xlim', [0 4]);
    hptheta=plot(iwdata.t, whisktheta, 'color', [0 175 255]/255, 'linewidth', 1);
    whiskrange=[min(whisktheta) max(whisktheta)];
    hold on
    if ~isempty(find(whisk_mask_ind))
        % make whisk_mask_ind patch
        % wh_episodes
        for iwh=1:size(wh_episodes, 1)
            p_wh(iwh) = patch(tvm([wh_episodes(iwh, 1) wh_episodes(iwh, 1) wh_episodes(iwh, 2) wh_episodes(iwh, 2)]), [whiskrange(1)  whiskrange(2) whiskrange(2)  whiskrange(1)],'r', 'facecolor', 'r','FaceAlpha',0.5, 'EdgeColor','none');
        end;
    end;
    
    if ~isempty(find(nw_mask_ind_new))
        nw_episodes=dissectmat(nw_mask_ind_new, 10);
        for inw=1:size(nw_episodes, 1)
            p_nw(inw) = patch(tvm([nw_episodes(inw, 1) nw_episodes(inw, 1) nw_episodes(inw, 2) nw_episodes(inw, 2)]), [whiskrange(1)  whiskrange(2) whiskrange(2)  whiskrange(1)],'m', 'facecolor', 'm','FaceAlpha',0.2, 'EdgeColor','none');
        end;
    end;
    
    if ~isempty(find(touch_mask_compute{iwdata.mainwid+1}))
        touch_episodes=dissectmat(touch_mask_compute{iwdata.mainwid+1}, 10);
        for it=1:size(touch_episodes, 1)
            p_touch(it) = patch(tvm([touch_episodes(it, 1) touch_episodes(it, 1) touch_episodes(it, 2) touch_episodes(it, 2)]), [whiskrange(1)  whiskrange(2) whiskrange(2)  whiskrange(1)], 'c', 'facecolor', 'c','FaceAlpha',0.2, 'EdgeColor','none');
        end;
    end;
    
    if ~isempty(licks)
        plot(licks, mean(get(gca, 'ylim')), 'go', 'markerfacecolor', 'g')
    end
    
    uistack(hptheta,'top')
    box off
    ylabel('Whisk pos (o)')
    title ([strrep(iwdata.cellname, '_', '-') ' trial#' num2str(itrialnum)])
    
    ha2=subplot(2, 1, 2);
    set(gca, 'nextplot', 'add', 'xlim', [0 4], 'ylim', vrange);
    
%     line([0 2.5], [spkthfixed spkthfixed],  'color', 'b');
    box off
    plot(tvm, iVmorg, 'k', 'linewidth', 1);
    hold on
    plot(tvm, iVmnoap, 'color', [.5 .5 .5], 'linewidth', 1);
    
    
    if ~isempty(find(whisk_mask_ind))
        % make whisk_mask_ind patch
        % wh_episodes
        for iwh=1:size(wh_episodes, 1)
            p_wh(iwh) = patch(tvm([wh_episodes(iwh, 1) wh_episodes(iwh, 1) wh_episodes(iwh, 2) wh_episodes(iwh, 2)]), [vrange(1)  vrange(2) vrange(2)  vrange(1)],'r', 'facecolor', 'r','FaceAlpha',0.5, 'EdgeColor','none');
        end;
    end;
    
    if ~isempty(find(nw_mask_ind_new))
        nw_episodes=dissectmat(nw_mask_ind_new, 10);
        for inw=1:size(nw_episodes, 1)
            p_nw(inw) = patch(tvm([nw_episodes(inw, 1) nw_episodes(inw, 1) nw_episodes(inw, 2) nw_episodes(inw, 2)]), [vrange(1)  vrange(2) vrange(2)  vrange(1)],'m', 'facecolor', 'm','FaceAlpha',0.2, 'EdgeColor','none');
        end;
    end;
    
    if ~isempty(find(touch_mask_compute{iwdata.mainwid+1}))
        touch_episodes=dissectmat(touch_mask_compute{iwdata.mainwid+1}, 10);
        for it=1:size(touch_episodes, 1)
            p_touch(it) = patch(tvm([touch_episodes(it, 1) touch_episodes(it, 1) touch_episodes(it, 2) touch_episodes(it, 2)]), [vrange(1)  vrange(2) vrange(2)  vrange(1)], 'c', 'facecolor', 'c','FaceAlpha',0.2, 'EdgeColor','none');
        end;
    end;
    
    plot(tvm, count_mask_ind*5+vrange(1), 'color', [.7 .7 .7], 'linewidth', 5); 
    
    ylabel('Vm (mV)')
    xlabel('Time (s)')
    linkaxes([ha1, ha2], 'x')
    % ----------------------------------------
    %%
    
    
    if ismac
        pathname= [finddropbox '/Work/Experiment_results/EpochsAnalysis/Traces/'];
    else
        pathname='C:\Work\Dropbox\Work\Experiment_results\EpochsAnalysis\Traces\';
        
    end;
    cd (pathname)
    
    if ~exist(['Epochs_' iwdata.cellname],'dir')
        mkdir(['Epochs_' iwdata.cellname])
    end;
    cd (['Epochs_' iwdata.cellname])
    print (hf(i), ['Epochs_', iwdata.cellname '_trial' num2str(itrialnum)], '-dpng')
    close(hf(i))
end;
end;



% Epochs contain all data regarding behavioral epochs, including phase and
% angle this can be used for future analysis.
% 
% % this is the PSP amplitude data:
% dVm.cell=iwdata.cellname;
% dVm.params=params;
% dVm.NW=PSP_nw_all;
% dVm.WH=PSP_wh_all;
% dVm.Touch=PSP_touch_all;
% dVm.mainwid=iwdata.mainwid; % 0 or 1,,
% 
% dtt.NW=dtt_nw_all;
% dtt.WH=dtt_wh_all;
% dtt.Touch=dtt_TOUCH_all;
% dtt.mainwid=iwdata.mainwid;

% transitions
transit.cell=iwdata.cellname;
transit.Vm=Vmtransit;
transit.Vm2=Vmtransit2;
transit.Spk=Spktransit;
transit.tvm=[-0.2*Fs : 0.5*Fs]/Fs;
transit.whiskpos=Whiskpostransit;
transit.whiskamp=Whisktransit;
transit.twhisk=[-0.2*1000 : 0.5*1000]/1000;
transit.transit_index=transit_index;

% whisking stuff
Epochs.cell=iwdata.cellname;

% dataout contains all
dataout.cell=iwdata.cellname;
% dataout.dVm=dVm;
% dataout.dtt=dtt;
dataout.transit=transit;
dataout.Epochs=Epochs;

path=finddropbox;
pathname=[path '/Work/Experiment_results/EpochsAnalysis/Results/'];


if isempty(params.laser)
    save ([pathname 'dataout_' strrep(iwdata.cellname, '-', '_') name  '.mat'], 'dataout');
else
    save ([pathname 'dataout_' strrep(iwdata.cellname, '-', '_') name 'stim' params.laser '.mat'], 'dataout');
end;

close all;
