function val = curveShorteningCoupledAttempt2Forcing2( curve, field, x )

deriv=periodicSplineCurveDerivEval(curve, x);
normed=norm(deriv);


val=-(2*perbspeval(field,x)-1)*(1/normed)*(deriv(1));

end


