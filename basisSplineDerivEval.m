function val = basisSplineDerivEval( U, u, i, p )
% Return the derivative of the ith degree-p basis spline at u
% INPUT:
% U=knot vector
% u=eval point or points
% i=i'th basis spline
% p=degree

val=zeros(numel(u), 1);
for j=1:numel(u)
    val(j)=(p/(U(i+p)-U(i)))*basisSplineEval(U, u(j), i, p-1) - (p/(U(i+p+1)-U(i+1)))*basisSplineEval(U, u(j), i+1, p-1);
end


end

