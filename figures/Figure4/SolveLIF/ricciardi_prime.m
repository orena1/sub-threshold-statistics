function rate_p=ricciardi_prime(mu,sigma,tau,VT,Vr,eps)

rate_p=(ricciardi(mu+eps,sigma,tau,VT,Vr)-ricciardi(mu,sigma,tau,VT,Vr))./eps;