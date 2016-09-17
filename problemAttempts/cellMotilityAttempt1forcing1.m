function vals = cellMotilityAttempt1forcing1( curve, field, k1, x )
% Forcing term f=k_1*a_1

deriv=periodicSplineCurveDerivEval(curve, x);
l = size(deriv);

norm = ones(1,l(2));
fieldEvals=perbspeval(field, x);
vals = zeros(1,l(2));

for i=1:l(2)
    norm(i)= sqrt(deriv(1,i)*deriv(1,i) + deriv(2,i)*deriv(2,i));
    vals(i)=-(k1*fieldEvals(i))*(1/norm(i))*(-deriv(2,i));
end

end

