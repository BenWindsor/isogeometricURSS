function val = periodicBasisEval( U, u, i, p )
% Returns the value of the i'th degree-p periodic basis spline on the knot vector U at
% point u
% INPUT:
% U = knot vector
% u = eval point
% i = ith periodic B-spline
% p = degrree

normalVals=localBasisSplineVectorEval(U, u, i, p);

%Assign correct element accounting for edge cases of repeated knots
%On the left
if i<=p
    elem=1;
%On the right
elseif i>=(numel(U)-p)
    elem=numel(U)-(p+1);
%Normal middle cases
else
    elem=i-p;
end

operator=localPeriodicOperator(U, u, p, elem);
newVals=operator*normalVals;

%IMPLEMENT: adjust which value is returned depending on which curve is
%needed
%First spline should be the one starting fully at first knot
if elem==1
    val=newVals(3);
elseif elem==2
    val=newVals(1);
else 
    val=newVals(3);
end
end

