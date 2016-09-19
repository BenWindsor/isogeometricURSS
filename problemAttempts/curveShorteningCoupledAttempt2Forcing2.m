function val = curveShorteningCoupledAttempt2Forcing2( curve, field, x )

deriv=periodicSplineCurveDerivEval(curve, x);
l = size(deriv);

norm = ones(1,l(2));
fieldEvals=perbspeval(field, x);
val = zeros(1,l(2));

for i=1:l(2)
    norm(i)= sqrt(deriv(1,i)*deriv(1,i) + deriv(2,i)*deriv(2,i));
    val(i)=-(2*fieldEvals(i)-1)*(1/norm(i))*(deriv(1,i));
end


end


