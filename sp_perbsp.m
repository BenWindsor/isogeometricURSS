function sp = sp_perbsp( perbsp, msh )
% An analogue of sp_nurbs for making a periodic basis spline space of
% functions.

transform = 'grad-preserving';
knots   = perbsp.knots;
degree  = perbsp.order - 1;

sp = sp_scalar_perbsp (knots, degree, msh, transform);

end

