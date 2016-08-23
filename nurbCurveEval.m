function val = nurbCurveEval( U, u, p, ctrl, W )
% Returns the value of the nurb curve at point u 
% INPUT:
% U = knot vector
% u = eval point/s
% p = degrree
% ctrl = control points
% W = weights

%OUTPUT:
% numerator = [x1 x2 ...;
%              y1 y2 ...;] etc.

%in nurbs book n=m-p-1 however there the knots run from 0,...,m but here
%they run from 1,...,m so we need to use n=(m-1)-p-1
n=numel(U)-1-p-1;

numerator=zeros(numel(ctrl(:,1)), 1);
denominator=weightFunction(U, u, p, W);

for i=1:n
    numerator(:,1) = numerator(:,1) + ctrl(:,i)*W(i)*basisSplineEval(U, u, i, p);
end

val=numerator/denominator;

end

