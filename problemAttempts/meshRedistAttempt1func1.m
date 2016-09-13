function val = meshRedistAttempt1func1( curve, x )

deriv=periodicSplineCurveDerivEval(curve, x);
normed=norm(deriv);

val=(1/(normed*normed))*(deriv(2)*deriv(2));

end

