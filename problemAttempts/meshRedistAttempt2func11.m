function val = meshRedistAttempt2func11( curve, x )

deriv=periodicSplineCurveDerivEval(curve, x);
normed=norm(deriv);

val=(1/(normed*normed))*(deriv(2)*deriv(2));

end

