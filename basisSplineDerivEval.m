function val = basisSplineDerivEval( U, u, i, p )
% Return the derivative of the ith degree-p basis spline at u
% INPUT:
% U=knot vector
% u=eval point
% i=i'th basis spline
% p=degree

val=(p/(U(i+p)-U(i)))*basisSplineEval(U, u, i, p-1) - (p/(U(i+p+1)-U(i+1)))*basisSplineEval(U, u, i+1, p-1);


end

