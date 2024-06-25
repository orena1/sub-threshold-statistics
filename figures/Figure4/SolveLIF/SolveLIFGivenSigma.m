function [mu_out, del_out]=SolveLIFGivenSigma(NetParm,StatVariables,sig,mu_0, del_0)
% g=.1; % mSim 0.1 * 10-3 Simens ;  Sim*Volt=Amper
% %  mV*mSim=microAmper
% %dv/dt=g/C*V
% % I nA = mV/msec
% C=1; %nFarad
% tau=C/g;% in mSec
% Vrest=-55;Theta=-30;Vr=-70;
% %Vrest=0;Theta=15;Vr=0;

g=NetParm.g;C=NetParm.C;
Vrest=NetParm.Vrest;Theta=NetParm.Theta;Vr=NetParm.Vr;
tau=C/g;
sigTheta=NetParm.sigTheta;
muTheta=NetParm.muTheta;
% mu_0=2.2; sig_0=1.02; del_0=0.02;
% mu_0=0.2; del_0=0.42;


barNu=StatVariables.barNu;barNu2=StatVariables.barNu2;

D=@(z)exp(-((z.^2)/2))./(sqrt(2*pi));
% r=I/J_0;

a0=[mu_0,del_0];

if sigTheta==0
    sol=fsolve(@F,a0);
else
    sol=fsolve(@F_Thresh,a0);
end
%sol=lsqnonlin(@F,a0);


% [a,fval,exitflag,output,jacobian] =fsolve(@F,a0);
mu_out=sol(1);del_out=sqrt(sol(2));



    function y=F(a)
        mu=a(1); Del2=a(2);
        Del=sqrt(Del2);
        
        
        Frz = @ (z)  phi_r(mu+Del*z,sig).*D(z);
        Fr2z = @ (z)  (phi_r(mu+Del*z,sig).^2).*D(z);
        mFR=integral(Frz,-Inf,Inf);
        mFR2=integral(Fr2z,-Inf,Inf);
        
        y(1)=mFR-barNu;
        y(2)=(mFR2-mFR.^2)-barNu2;
        
        
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

    function y=F_Thresh(a)
        mu=a(1); Del2=a(2);
        
        Del=sqrt(Del2);
        
        Frz = @ (z,x)  phi_r(mu+Del*z,sig,x).*D(z);
        Fr2z = @ (z,x)  (phi_r(mu+Del*z,sig,x).^2).*D(z);
        % for uniform distribution 1/(2a) int_-a^a dx whre a=sqrt(3)*sd
        mFR=integral2(Frz,-Inf,Inf,muTheta-sqrt(3)*sigTheta,muTheta+sqrt(3)*sigTheta)./(2*sqrt(3)*sigTheta);
        mFR2=integral2(Fr2z,-Inf,Inf,muTheta-sqrt(3)*sigTheta,muTheta+sqrt(3)*sigTheta)./(2*sqrt(3)*sigTheta);
        
        y(1)=mFR-barNu;
        y(2)=(mFR2-mFR.^2)-barNu2;
        
        function Fr=phi_r(I,sig,thet)
            Fr=ricciardi((I)/g+Vrest,sig/sqrt(g*C),tau,thet,Vr);
        end
        
        function v=phi_v(I,sig,Theta)
            v=Vrest+(I)./g-(thet-Vr).*tau.*phi_r(I,sig,thet);
        end
        
    end


end