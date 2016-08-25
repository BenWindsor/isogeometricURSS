function val = periodicBasisEval( U, u, i, p )
% Returns the value of the i'th degree-p periodic basis spline on the knot vector U at
% point u
% INPUT:
% U = knot vector
% u = eval point
% i = ith periodic B-spline
% p = degrree

normalVals=localBasisSplineVectorEval(U, u, i, p);
elem=i-p;
operator=localPeriodicOperator(U, u, p, elem);
newVals=operator*normalVals;
val=newVals(2);

end

