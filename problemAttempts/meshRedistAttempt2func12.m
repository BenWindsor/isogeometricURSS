function val = meshRedistAttempt2func12( curve, x )

deriv=periodicSplineCurveDerivEval(curve, x);
l = size(deriv);
normsqr = ones(1,l(2));
for j=1:1:l(2)
    normsqr(j) = deriv(:,j).*deriv(:,j);
end

val=(-1/(normed*normed))*(deriv(1)*deriv(2));

end

