function val = meshRedistAttempt2func22( curve, x )

deriv=periodicSplineCurveDerivEval(curve, x);
normed=norm(deriv);

val=(1/(normed*normed))*(deriv(1)*deriv(1));

end

