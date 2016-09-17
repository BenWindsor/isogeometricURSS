function vals = cellMotilityAttempt1lambdaFunc1( curve, x )

deriv=periodicSplineCurveDerivEval(curve, x);
l = size(deriv);

norm = ones(1,l(2));
vals = zeros(1,l(2));

for i=1:l(2)
    norm(i)= sqrt(deriv(1,i)*deriv(1,i) + deriv(2,i)*deriv(2,i));
    vals(i)= (1/norm(i))*(-deriv(2,i));
end

end

