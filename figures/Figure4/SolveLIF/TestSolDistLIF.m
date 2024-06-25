function [barNu, barNu2, barV,barV2, barSigV, barSigV2]=TestSolDistLIF(mu, sig, Del,NetParm)

% g=0.1;C=1;tau=C/g;
% %Vrest=0;Theta=15;Vr=0;
%  Vrest=-55;Theta=-30;Vr=-70;

g=NetParm.g;C=NetParm.C;
Vrest=NetParm.Vrest;Theta=NetParm.Theta;Vr=NetParm.Vr;
tau=C/g;

D=@(z)exp(-((z.^2)/2))./(sqrt(2*pi));

y=F(mu,sig,Del);
barNu=y(1);barNu2=y(2);barV=y(3);barV2=y(4);
barSigV=y(5); barSigV2=y(6);

% barNu=100*barNu;
% barNu2=100^2*barNu2;
    function y=F(mu,sig,Del)
        
        Frz = @ (z)  phi_r(mu+Del*z,sig).*D(z);
        Fr2z = @ (z)  (phi_r(mu+Del*z,sig).^2).*D(z);
        
        Vz = @ (z)  phi_v(mu+Del*z,sig).*D(z);
        V2z = @ (z)  (phi_v(mu+Del*z,sig).^2).*D(z);
        
        mFR=integral(Frz,-Inf,Inf);
        mFR2=integral(Fr2z,-Inf,Inf);
        
        mV=integral(Vz,-Inf,Inf);
        mV2=integral(V2z,-Inf,Inf);
        
        y(1)=mFR;
        y(2)=(mFR2-mFR.^2);
        y(3)=mV;
        y(4)=(mV2-mV.^2);
        
        sigVz = @ (z)  phi_sigV(mu+Del*z,sig).*D(z);
        sigV2z = @ (z)  (phi_sigV(mu+Del*z,sig).^2).*D(z);
        
        msigV=integral(sigVz,-Inf,Inf);
        msigV2=integral(sigV2z,-Inf,Inf);
        
        y(5)=msigV;
        y(6)=(msigV2-msigV.^2);  
    end


     function Fr=phi_r(I,sig)
        Fr=ricciardi((I)/g+Vrest,sig/sqrt(g*C),tau,Theta,Vr);
    end

    function v=phi_v(I,sig)
        v=Vrest+(I)./g-(Theta-Vr)*tau*phi_r(I,sig);
    end


    function sigV=phi_sigV(I,sig)
        sigV=(tau/2*(sig/C)^2+(Theta-Vr).*(phi_v(I,sig)-(Theta+Vr)./2 ).*tau.*phi_r(I,sig));
    end


end