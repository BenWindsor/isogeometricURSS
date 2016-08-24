function val = nurbBasisEval( U, u, i, p, W )
% Returns value of nurbs basis function R_i,p(u) 
% INPUTS
% U=knot vector
% u=eval point
% i=the ith nurbs basis
% p=degree
% W=weight vector

numerator=basisSplineEval(U,u,i,p)*W(i);
denominator=weightFunction(U,u,p,W);

val=numerator/denominator;

end

