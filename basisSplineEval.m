function val = basisSplineEval( U, u, i, p )
% Returns the value of the i'th degree p spline on the knot vector U at
% point u
% INPUT:
% U = knot vector
% u = eval point
% i = ith B-spline
% p = degrree

if p==0
    val = zeroDegreeEval(U, u, i);
else 
    val = firstCoeff(U, u, i, p)*basisSplineEval(U, u, i, p-1) + secondCoeff(U, u, i, p)*basisSplineEval(U, u, i+1, p-1);
end
    