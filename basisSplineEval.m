function val = basisSplineEval( U, u, i, p )
% Returns the value of the i'th degree p spline on the knot vector U at
% point u
% INPUT:
% U = knot vector
% u = eval point/s
% i = ith B-spline
% p = degrree

val=zeros(numel(u),1);

for j=1:numel(u)
    if p==0
        val(j) = zeroDegreeEval(U, u(j), i);
    else
        val(j) = firstCoeff(U, u(j), i, p)*basisSplineEval(U, u(j), i, p-1) + secondCoeff(U, u(j), i, p)*basisSplineEval(U, u(j), i+1, p-1);
    end
end
end
    