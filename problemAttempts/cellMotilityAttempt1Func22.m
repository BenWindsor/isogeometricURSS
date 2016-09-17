function vals = cellMotilityAttempt1Func22( curve, x)

deriv=periodicSplineCurveDerivEval(curve, x);
l = size(deriv);

normsqr = ones(1,l(2));
vals = zeros(1,l(2));

for j=1:1:l(2)
    normsqr(j) = deriv(1,j)*deriv(1,j) + deriv(2,j)*deriv(2,j);
    vals(j) = (1/normsqr(j))*(deriv(1,j)*deriv(1,j));
end

end

