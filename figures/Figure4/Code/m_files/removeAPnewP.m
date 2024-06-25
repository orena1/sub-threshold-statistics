function [Vm_removed, spikes]=removeAPnewP(Vm,Fs,ratio_dvdt, Vbound, dvdtlow, toplot)
% removeAPnew(Vm,Fs,0.33,  [-50 -40 50 10], 10000)
% Vth=detectSpikeonsetVm(Vm, 0.33, [-50 -40 50 10], 10000, 10000);
%  ratio_dvdt, vbound, dvdtlow,
% function Vm_removed=removeAP(Vm,Fs,threshold, width)
% 9.11.2014 revised to deal with bursts (inter-spike interval<10 ms, instantaneous rate>100 Hz

% vth.Vboundtype
%     'threhold_lowest_possible'
%     'peak_lowest_possible'
%     'peak_highest_possible'
%     'size_smallest_possible'
%     e.g., vbound=[[-50 -40 50 10]]
%     vth=detectSpikeonsetVm(vx, 0.33, [-50 -40 50 10], 10000, 10000);

if nargin<6
    toplot=0;
end;

if isempty(Vbound)
    Vbound=[-50 -40 50 10];
end;

% detectSpikeonsetParray(P, ratio_dvdt, vbound, dvdtlow,  checkAP, toplot)
Vth=detectSpikeonsetVmP(Vm, ratio_dvdt, Vbound, dvdtlow, Fs, toplot);
spikes=Vth.spktrain;
[L, repeats]=size(Vm);

% pre_spike=round(0.5*Fs/1000);  % number of data points before spike threshold
% post_spike=round(4*Fs/1000);  % number of data points after spike threshold

pre_spike=round(20*Fs/1000);  % number of data points before spike threshold
post_spike=round(20*Fs/1000);  % number of data points after spike threshold

% Changed on May 29 2018
pre_spike=round(5*Fs/1000);  % number of data points before spike threshold
post_spike=round(5*Fs/1000);  % number of data points after spike threshold

pre_spike=round(1*Fs/1000);  % number of data points before spike threshold
post_spike=round(1*Fs/1000);  % number of data points after spike threshold

V_mask=zeros(size(Vm));
Vm_removed=zeros(size(Vm));
for itrial=1:repeats
    timings=find(Vth.spktrain(:,itrial));
    if any(timings)  % means there are spikes
        V_unremoved=Vm(:,itrial);
        
        for ispk=1:length(timings)
            V_mask(max([1 timings(ispk)-pre_spike]):min([timings(ispk)+post_spike, L]), itrial)=1;
        end;
        % if the gap between "clipps" is less than 2 ms, linked them
        % [eps, newdatamat]=dissectmat(datamat, mindur)
        [V_mask(:, itrial)]=mergegaps(V_mask(:, itrial), 2*Fs/1000, 2);
        eps=dissectmat(V_mask(:, itrial), 1);
        % eps is the begs and ends of each region that needs to be clipped.
        Vm_removed(V_mask(:, itrial)==0, itrial)=Vm(V_mask(:, itrial)==0,itrial);
        % then fill the rest
        
        for k=1:size(eps)
            
            if eps(k, 1)==1 % odd ball, start at the very beginning
                Vm_removed(eps(k, 1):eps(k, 2), itrial)=Vm(eps(k, 2), itrial);
            elseif eps(k, 2)==L % start at the very end
                Vm_removed(eps(k, 1):eps(k, 2), itrial)=Vm(eps(k, 1), itrial);
            else % most common situation
                Vm_removed(eps(k, 1):eps(k, 2), itrial)=interp1([eps(k, 1) eps(k, 2)], Vm([eps(k, 1) eps(k, 2)], itrial), [eps(k, 1):eps(k, 2)]);
                % interp1([n_pre n_post],[V_pre' V_post'],[n_pre:n_post]);
            end
            
        end;
        
        
    else
        Vm_removed(:, itrial)=Vm(:, itrial);
    end;
    if toplot
        figure(34); clf
        plot(Vm(:, itrial), 'k'); hold on
        plot(Vm_removed(:, itrial), 'b');
        pause (1)
    end;
end;