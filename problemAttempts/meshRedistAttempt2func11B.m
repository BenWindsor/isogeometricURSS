function val = meshRedistAttempt2func11B( curve, x )

deriv=periodicSplineCurveDerivEval(curve, x);
l = size(deriv);

normsqr = ones(1,l(2));
val = zeros(1,l(2));

for j=1:1:l(2)
    normsqr(j) = deriv(1,j)*deriv(1,j) + deriv(2,j)*deriv(2,j);
    val(j) = (1/normsqr(j))*(deriv(2,j)*deriv(2,j));
end

%val=val';
end



