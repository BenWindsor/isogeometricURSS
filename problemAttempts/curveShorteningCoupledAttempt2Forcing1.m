function val = curveShorteningCoupledAttempt2Forcing1( curve, field, x )

deriv=periodicSplineCurveDerivEval(curve, x);
normed=norm(deriv);

val=-(2*perbspeval(field, x)-1)*(1/normed)*(-deriv(2));
end

