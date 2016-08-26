function varargout = geo_perbspline (perbspline, deriv, deriv2, pts, ders, rdim)

% construct a geometry map from the periodic basis spline structure
%
%   output = geo_nurbs (perbspline, dnurbs, dnurbs2, pts, ders)
%   output = geo_nurbs (perbspline, dnurbs, dnurbs2, pts, ders, rdim)
%
% INPUTS:
%
%   nurbs:    periodic basis spline structure that defines the geometry
%   dnurbs:   periodic basis spline  structure for the derivatives (see pernrbderiv.m)
%   dnurbs2:  periodic basis spline  structure for the second derivatives (see perbspderiv.m)
%   pts    :  points where the map has to be evaluated
%   ders   :  number of derivatives to be evaluated (from 0 to 2)
%   rdim   :  the dimension of the physical domain
%   
% OUTPUT:
%
%   output: for ders = 0, the parametrization F evaluated at pts
%           for ders = 1, the Jacobian of the parametrization, evaluated at pts
%           for ders = 2, the Hessian of the parametrization, evaluated at pts
%
%   Multiple outputs are also allowed, for more efficient computations. For instance
%
%     [F, jac] = geo_perbspline (perbspline, pts, 1);
%     [F, jac, hess] = geo_perbspline (perbspline, pts, 2);


  ndim = numel (perbspline.order);

  if (nargin < 4)
    if (any (abs(perbspline.coefs(3,:)) > 1e-12))
      rdim = 3;
    elseif (any (abs(perbspline.coefs(2,:)) > 1e-12))
      rdim = 2;
    else
      rdim = 1;
    end
  end

  if (rdim < ndim)
    error ('geo_perbspline: the dimensions of the geometry are wrong. rdim < ndim?')
  end
  
  switch (ders)
    case 0
      F = nrbeval (perbspline, pts);
      varargout{1} = F(1:rdim, :);
    case 1
%       deriv = nrbderiv (nurbs);
      [F, jac] = nrbdeval (perbspline, deriv, pts);
      if (iscell (pts))
        npts = prod (cellfun (@numel, pts));
        map_jac = zeros (rdim, ndim, npts);
        for idim = 1:ndim
          map_jac(1:rdim, idim, :) = reshape (jac{idim}(1:rdim,:,:), rdim, 1, npts);
        end
      else
        map_jac = zeros (rdim, ndim, size (pts, 2));
        for idim = 1:ndim
          map_jac(1:rdim, idim, :) = reshape (jac{idim}(1:rdim, :), rdim, 1, size (pts, 2));
        end
      end
      if (nargout == 1)
        varargout{1} = map_jac;
      elseif (nargout == 2)
        varargout{1} = F(1:rdim, :);
        varargout{2} = map_jac;
      end

    case 2
%       [deriv, deriv2] = nrbderiv (nurbs);
      [F, jac, hessian] = nrbdeval (perbspline, deriv, deriv2, pts);

      if (iscell (pts))
        npts = prod (cellfun (@numel, pts));
        map_jac = zeros (rdim, ndim, npts);
        hess = zeros (rdim, ndim, ndim, npts);
        for idim = 1:ndim
          map_jac(1:rdim, idim, :) = reshape (jac{idim}(1:rdim,:,:), rdim, 1, npts);
          for jdim = 1:ndim
            hess(1:rdim, idim, jdim, :) = reshape (hessian{idim,jdim}(1:rdim,:,:), rdim, 1, 1, npts);
          end
        end
      else
        map_jac = zeros (rdim, ndim, size (pts, 2));
        hess = zeros (rdim, ndim, ndim, size (pts, 2));
        for idim = 1:ndim
          map_jac(1:rdim, idim, :) = reshape (jac{idim}(1:rdim, :), rdim, 1, size (pts, 2));
          for jdim = 1:ndim
            hess(1:rdim, idim, jdim, :) = reshape (hessian{idim,jdim}(1:rdim,:), rdim, 1, 1, size (pts, 2));
          end
        end
      end
      if (nargout == 1)
        varargout{1} = hess;
      elseif (nargout == 3)
        varargout{1} = F(1:rdim, :);
        varargout{2} = map_jac;
        varargout{3} = hess;
      end
    otherwise
      error ('geo_nurbs: number of derivatives limited to two')
  end

end
