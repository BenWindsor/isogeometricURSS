function val = firstCoeff( U, u, i, p )
% Returns the value of the first coefficient in the spline recursive
% definition
% INPUT:
% U = knot vector
% u = eval point
% i = ith B-spline
% p = degrree

if (U(i+p)-U(i))==0
    val = 0;
else
    val=(u-U(i))/(U(i+p)-U(i));
end

end

