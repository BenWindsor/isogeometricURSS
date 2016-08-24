function val = weightFunction( U, u, p, W )
% Returns the value of the NURBS weight function at u
% U = knot vector
% u = eval point/s
% p = degrree
% W = weight vector

% This function is the denominator in the NURBS basis functions

n = numel(U)-p-1;
val = 0;
for i=1:n
    val = val + W(i)*basisSplineEval(U, u, i, p);
end

end

