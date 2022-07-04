COMMENT
Two state kinetic scheme synapse described by rise time tau1,
and decay time constant tau2. The normalized peak conductance is 1.
Decay time MUST be greater than rise time.

The solution of A->G->bath with rate constants 1/tau1 and 1/tau2 is
 A = a*exp(-t/tau1) and
 G = a*tau2/(tau2-tau1)*(-exp(-t/tau1) + exp(-t/tau2))
	where tau1 < tau2

If tau2-tau1 is very small compared to tau1, this is an alphasynapse with time constant tau2.
If tau1/tau2 is very small, this is single exponential decay with time constant tau2.

The factor is evaluated in the initial block 
such that an event of weight 1 generates a
peak conductance of 1.

Because the solution is a sum of exponentials, the
coupled equations can be solved as a pair of independent equations
by the more efficient cnexp method.

ENDCOMMENT

NEURON {
	POINT_PROCESS Exp2Syn_cc
	RANGE tau1, tau2, e, i, alpha, rho
	NONSPECIFIC_CURRENT i

	RANGE g, e_pas
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau1 = 0.1 (ms) <1e-9,1e9>
	tau2 = 10 (ms) <1e-9,1e9>
	e=0	(mV)
	alpha =  1
	rho = 1
	e_pas = -900000 (mV) :default value is bad!
}

ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
	factor
}

STATE {
	A (uS)
	B (uS)
}

INITIAL {
	LOCAL tp
	if (tau1/tau2 > 0.9999) {
		tau1 = 0.9999*tau2
	}
	if (tau1/tau2 < 1e-9) {
		tau1 = tau2*1e-9
	}
	A = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = B - A
	i = g*alpha*(rho*(v - e)+(1-rho)*(e_pas-e))
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}

NET_RECEIVE(weight (uS)) {
	A = A + weight*factor
	B = B + weight*factor
}
