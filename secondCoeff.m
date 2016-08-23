function val = secondCoeff( U, u, i, p )
% Returns the value of the second coefficient in the spline recursive
% definition
% INPUT:
% U = knot vector
% u = eval point
% i = ith B-spline
% p = degrree

if (U(i+p+1)-U(i+1))==0
    val = 0;
else
    val=(U(i+p+1)-u)/(U(i+p+1)-U(i+1));
end

end

