function val = periodicSingleBasisDerivEval( U, u, elem, p )
% Takes a single periodic basis and returns its derivative at a point.
% An analogue of the bspderiv.m file in the nurbs toolbox for periodic
% basis functions
% INPUT
% U=UNIFORM OPEN knot vector
% u=eval point
% elem=element the periodic basis spline starts at
% p=degree

% Method: as all periodic basis splines are the same, they match with one
% of the normal B-splines just shifted. We compute the derivative at a
% point on the B-Spline that corresponds to our point on the periodic basis
% spline.

if (p>3)
    error('Not implemented for p>3');
end

index=elem+p; %index in the knots where that element begins
offset=periodicOffset(U, U(index), u); %How far from the start of that knot u is periodically
val=basisSplineDerivEval(U, U(p+1)+offset, p+1, p); %number p+1 is first spline to be like periodic splines

end

