function val = periodicSplineSurfEval( U, u, p, ctrl )
% Returns the value of the spline curve at point u 
% INPUT:
% U = OPEN UNIFORM knot vector
% u = cell of eval point/s
% p = degrree HAS TO BE DEGREE 2 or 3 at the moment
% ctrl = cell of control points in dimensions higher then structure ctrl like:
% [C1,1 C1,2 C1,3 C1,4;
%  C2,1 C2,2 C2,3 C2,4] etc.
% So each row is a row of control points in one direction 

val=cell(numel(u),1);
elems=numel(U)-(2*p)-1;

for ipoint=1:numel(u)
    point=u{ipoint};
    val{ipoint}=[0];
    for i=1:elems
        for j=1:elems
            val{ipoint}=val{ipoint}+periodicBasisEval(U, point(1), i, p)*periodicBasisEval(U, point(2), i, p)*ctrl(i,j);
        end
    end
end

end


