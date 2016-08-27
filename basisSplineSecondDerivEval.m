function val = basisSplineSecondDerivEval( U, u, i, p)
% Return the second derivative of the ith degree-p basis spline at u
% INPUT:
% U=knot vector
% u=eval point or points
% i=i'th basis spline
% p=degree

%Using 'The NURBS Book' algo pg61

val=zeros(numel(u), 1);
for j=1:numel(u)
    val(j)=p*((basisSplineDerivEval(U, u(j), i, p-1)/(U(i+p)-U(i))) - (basisSplineDerivEval(U, u(j), i+1, p-1)/(U(i+p+1)-U(i+1))));
end

end

