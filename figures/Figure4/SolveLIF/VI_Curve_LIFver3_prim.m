function [v, sigV2, Fr,Fr_prim]=VI_Curve_LIFver3_prim(mu,sig,Del,N,NetParm)

g=NetParm.g;C=NetParm.C;
Vrest=NetParm.Vrest;Theta=NetParm.Theta;Vr=NetParm.Vr;
tau=C/g;
sigTheta= NetParm.sigTheta;
muTheta=    NetParm.muTheta;

%  g=.1; % mSim 0.1 * 10-3 Simens ;  Sim*Volt=Amper
%  %  mV*mSim=microAmper
%  %dv/dt=g/C*V
%  % I nA = mV/msec
%  C=1; %nFarad
%  tau=C/g;% in mSec
% % % % %   Vrest=0;Theta=15;Vr=0;
%   Vrest=-55;Theta=-30;Vr=-70;



% g=1;C=1;tau=C/g;tau=0.1;
% Vrest=0;Theta=1;Vr=0;

%Vrest=-65;Theta=-30;Vr=-70;

I=mu+Del*randn(N,1);
if sigTheta==0
    
    [v, sigV2,  Fr,Fr_prim]=phiV(I,sig);
    
    
else
    thet=muTheta-sqrt(3)*sigTheta+2*sqrt(3)*sigTheta.*(rand(N,1));
    [v, sigV2,  Fr]=phiVThresh(I,sig,thet);
end





    function [v,sigV2, Fr,Fr_prim]=phiV(I,sig)
        %          Xmax=sqrt(2*g*C)/sig*(Theta-((I)/g+Vrest));
        %          Xmin=sqrt(2*g*C)/sig*(Vr-((I)/g+Vrest));
        %         if Xmax>20
        %              y=Vrest+(I)./g-(Theta-Vr)*tau*g/C*(Xmax/sqrt(2*pi)*exp(-Xmax^2));
        %         else
        %         y=Vrest+(I)./g-(Theta-Vr)*tau*g/C*(sqrt(2*pi)*F(Xmin,Xmax))^(-1);
        %         end
        
        
        Fr=ricciardi((I)/g+Vrest,sig/sqrt(g*C),tau,Theta,Vr);
        Fr_prim=ricciardi_prime((I)/g+Vrest,sig/sqrt(g*C),tau,Theta,Vr,0.01);
        %           Fr=ricciardi((I)/g+Vrest,sig/sqrt(2*g*C),tau,Theta,Vr);% WHY 100????why 2
        
        %         if Xmax>20
        %              y=Vrest+(I)./g-(Theta-Vr)*(Xmax/sqrt(2*pi)*exp(-Xmax^2));
        %         else
        v=Vrest+(I)./g-(Theta-Vr)*tau*Fr;
        sigV2=tau/2*(sig/C)^2+(Theta-Vr).*(v-(Theta+Vr)./2 ).*tau.*Fr;
        %         end
        %         Fr=1/tau*(sqrt(2*pi)*FHV(Xmin,Xmax))^(-1);
        % Fr=1;
        
    end

    function [v,sigV2, Fr]=phiVThresh(I,sig,thet)
        Fr=ricciardi((I)/g+Vrest,sig/sqrt(g*C),tau,thet,Vr);
        v=Vrest+(I)./g-(thet-Vr).*tau.*Fr;
        sigV2=tau/2*(sig/C)^2+(thet-Vr).*(v-(thet+Vr)./2 ).*tau.*Fr;
        
    end



    function y=FHV(Xmin,Xmax)
        
        f=@ (x) exp(x.^2/2).*erfc(x);
        y=integral(f,Xmin,Xmax);
    end


%     function y=F(Xmin,Xmax)
%
%         f=@ (x) exp(x.^2/2).*erfc(-x);
%         y=integral(f,Xmin,Xmax);
%
%     end


end