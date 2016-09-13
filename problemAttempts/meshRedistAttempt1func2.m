function val = meshRedistAttempt1func2( curve, x )

deriv=periodicSplineCurveDerivEval(curve, x);
normed=norm(deriv);

val=(1/(normed*normed))*(deriv(1)*deriv(1));

end

