function [Tout,Yout,sp]=Test_LIF_DiffAprrox_Ver2(mu,sig,NetParm,SimParam)
% Bug 13/4/2020:  for c!=1 it genetrates a bug.
%    dv = (-(v-Vrest))./tau+I;   instead of      dv = (-(v-Vrest))./tau which is wrong;    
g=NetParm.g;C=NetParm.C;
Vrest=NetParm.Vrest;Vr=NetParm.Vr;
vthresh=NetParm.Theta;
N=NetParm.N;

dt=SimParam.dt;
t_end=SimParam.t_end;
t_0=SimParam.t_0;


tau=C/g;
I=mu./C; %Here there was a bug of mu/g instead of mu/C
% sig=sig./(g*tau);
sig=sig./(C);
Y0=Vr+randn(N,1);

% [Tout,Yout]= ode45(@fun_LIF,tspan,u0');
% [Tout,Yout]=euler(@fun_WB,tspan(1),tspan(end),u0',0.05);
[Tout,Yout,sp]= Euler_Maruyama(@fun_LIF,1,t_0,t_end,Y0,dt);


    function dv=fun_LIF(v)        
        dv = (-(v-Vrest))./tau+I;   
    end

    function [T,Y,sp]=Euler_Maruyama(f,g,t_0,t_end,Y0,dt)
        
        % Input:
        %        Set discretization [t_0:dt:t_end]
        %        Set initial condition Y0
        %        Set the functions f,g: dx=f(x)dt+g(x)dW
        % Output:
        %        T- time vector
        %        Y- estimate of Y(t)
        
        n=round((t_end-t_0)/dt);
        Y=zeros(length(Y0),n+1);
        T=t_0:dt:t_end;
        Y(:,1)=Y0;
        sp=[];
        
        dW=sig.*sqrt(dt)*randn(length(Y0),n); % Array of Brownian movements
        
        
        for i=1:n
            Y(:,i+1) = Y(:,i) + dt*f(Y(:,i)) + dW(:,i);
            if sum((Y(:,i+1))>vthresh)>0
                indSp=  (Y(:,i+1)>vthresh)   ;
                numSpikes=sum(indSp);
                sp =[sp; T(i).*ones(numSpikes,1),find(Y(:,i+1)>vthresh)];
                Y(indSp,i+1)=Vr(indSp);
            end
        end
    end

end


