function val = periodicSplineCurveEval( U, u, p, ctrl )
% Returns the value of the spline curve at point u 
% INPUT:
% U = OPEN UNIFORM knot vector
% u = eval point/s
% p = degrree HAS TO BE DEGREE 2 at the moment
% ctrl = control points in dimensions higher then structure ctrl like:
% [x1 x2 x3 x4;
%  y1 y2 y3 y4] etc.
% so each column is a coordinate point

%convert cell of points to normal array
if iscell(u)
    u=cell2mat(u);
end

val=zeros(numel(ctrl(:,1)), numel(u));
elems=numel(U)-2*p-1;
for i=1:numel(u)
    for j=1:elems
        val(:,i) = val(:,i) + ctrl(:,j)*periodicBasisEval(U, u(i), j, p);
    end
end

end
