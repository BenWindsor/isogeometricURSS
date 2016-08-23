function val = nurbBasisEval( U, u, i, p, W )
%NURBBASISEVAL Summary of this function goes here
%   Detailed explanation goes here

numerator=basisSplineEval(U,u,i,p)*W(i);
denominator=weightFunction(U,u,p,W);

val=numerator/denominator;

end

