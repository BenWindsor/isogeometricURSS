function vals = bernsteinVectorEval(u, p, n, m )
% Returns the vector of bernstein polynomials evaluated at the point u and
% of length n+m as in the bezier extraction paper just above eqn(18)
% INPUT
% u=eval point
% p=degree
% n=number of control points
% m=total number of knots required to produce bezier decomposition

% note: can obtain m from numel(new) in return values of bezierKnots.m
vals=zeros(n+m,1);

for i=1:(n+m)
    vals(i)=bernsteinEval(i,p,u);
end

end


