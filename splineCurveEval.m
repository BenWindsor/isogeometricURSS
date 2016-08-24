function val = splineCurveEval( U, u, p, ctrl )
% Returns the value of the spline curve at point u 
% INPUT:
% U = knot vector
% u = eval point/s
% p = degrree
% ctrl = control points in dimensions higher then structure ctrl like:
% [x1 x2 x3 x4;
%  y1 y2 y3 y4] etc.
% so each column is a coordinate point

val=zeros(numel(ctrl(:,1)), numel(u));
n=numel(U)-p-1
for i=1:numel(u)
    for j=1:n
        val(:,i) = val(:,i) + ctrl(:,j)*basisSplineEval(U, u(i), j, p)
    end
end

end

