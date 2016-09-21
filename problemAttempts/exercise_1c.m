% EXERCISE_1C: solution to exercise 1 c) of the tutorial, based on ex_article_15lines.
%
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011, 2012 Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Construct the geometry of exercise 1, a)
line = nrbline ([1 0], [2 0]);
srf = nrbrevolve (line, [0 0 0], [0 0 1], pi/2);
srf = nrbdegelev (srf, [0 1]);

% Set the number of knots that have to be inserted in the original geometry
%  to obtain the desired refined meshes. Remeber that, in exercise 1, b), 
%  for a 10x10 mesh we inserted 9 knots in each direction.
nsub = [0 1 3 7 15 31 63];

% Loop in the meshes
for isub = 1:numel (nsub)
% Compute the knots to be inserted in the coarse geometry. 
%  Degree 2 and continuity C^1, in both directions.
  [~, ~, nknots] = kntrefine (srf.knots, ...
                            [nsub(isub) nsub(isub)], [2 2], [1 1]);
% This is the alternative for C^0 continuity
% UNCOMMENT THE FOLLOWING LINES TO SOLVE WITH C^0 CONTINUITY
%  [~, ~, nknots] = kntrefine (srf.knots, ...
%                            [nsub(isub) nsub(isub)], [2 2], [0 0]);

% Perform the knot insertion, and save the result in a new NURBS structure.
  nurbs = nrbkntins (srf, nknots);
% The function geo_load also works with a NURBS structure (see the help).
  geometry = geo_load (nurbs);
% This part remains as in the original example.
  knots = geometry.nurbs.knots;

  [qn, qw] = msh_set_quad_nodes (knots, msh_gauss_nodes (geometry.nurbs.order));
  msh = msh_2d (knots, qn, qw, geometry);

  space  = sp_nurbs_2d (geometry.nurbs, msh);

  mat = op_gradu_gradv_tp (space, space, msh, @(x, y) ones (size (x)));
  rhs = op_f_v_tp (space, msh, @(x, y) (8-9*sqrt(x.^2+y.^2)).*sin(2*atan(y./x))./(x.^2+y.^2));

  drchlt_dofs = [];
  for iside = 1:4
    drchlt_dofs = union (drchlt_dofs, space.boundary(iside).dofs);
  end
  int_dofs = setdiff (1:space.ndof, drchlt_dofs);

  u = zeros (space.ndof, 1);
  u(int_dofs) = mat(int_dofs, int_dofs) \ rhs(int_dofs);

  err = sp_l2_error (space, msh, u, @(x,y)(x.^2+y.^2-3*sqrt(x.^2+y.^2)+2).*sin(2.*atan(y./x)));

% Save the number of degrees of freedom and the L2 error in two variables.
  ndof(isub) = space.ndof;
  err_l2(isub) = err
end

% Plot the error vs. the degrees of freedom in logarithmic scale.
loglog (ndof, err_l2, 'b-o', ndof, ndof.^(-1.5), 'k')