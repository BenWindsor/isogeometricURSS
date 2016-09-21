% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data
% The domain is the extrusion of geometry of exercise 1a).
line = nrbline ([1 0], [2 0]);
srf = nrbrevolve (line, [0 0 0], [0 0 1], pi/2);
srf = nrbdegelev (srf, [0 1]);
vol = nrbextrude (srf, [0 0 5]);
% We save the geometry in a txt file, and set it as the geometry of the problem.
nrbexport (vol, 'cylinder.txt')
problem_data.geo_name = 'cylinder.txt';

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [4 6];
problem_data.drchlt_sides = [];
problem_data.press_sides  = [3];
problem_data.symm_sides   = [1 2 5];

% Physical parameters, as in the 2D example (ex_plane_strain_ring), but the
%  function handles must also take the z coordinate as an argument.
E  =  1; nu = 0;
problem_data.lambda_lame = @(x, y, z) ((nu*E)/((1+nu)*(1-2*nu)) * ones (size (x)));
problem_data.mu_lame = @(x, y, z) (E/(2*(1+nu)) * ones (size (x)));

% Source and boundary terms. They are all zero, except the internal pressure.
%  They are as in the 2D example, but the z coordinate must be an argument.
P = 1;
problem_data.f = @(x, y, z) zeros (3, size (x, 1), size (x, 2));
problem_data.g = @(x, y, z, ind) zeros (3, size (x, 1), size (x, 2));
problem_data.h = @(x, y, z, ind) zeros (3, size (x, 1), size (x, 2));
problem_data.p = @(x, y, z, ind) P * ones (size (x));

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
%  We set the degree 3 and regularity C2 in the three directions,
%  5 elements in each direction (125 elements), and 4 quadrature points 
%  per element in each direction (64 points per element)
clear method_data
method_data.degree     = [3 3 3];     % Degree of the bsplines
method_data.regularity = [2 2 2];     % Regularity of the splines
method_data.nsub       = [5 5 5];     % Number of subdivisions
method_data.nquad      = [4 4 4];     % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
%  The solver is the one for 3D linear elasticity, as in 
%   the example ex_lin_elast_horseshoe.
[geometry, msh, space, u] = solve_linear_elasticity_3d (problem_data, method_data);

% 4) POST-PROCESSING. EXPORT TO PARAVIEW
%  The postprocessing is as in the 3D example (ex_lin_elast_horseshoe).
output_file = 'lin_elast_cylinder_Deg3_Reg2_Sub1';

vtk_pts = {linspace(0, 1, 30), linspace(0, 1, 30), linspace(0, 1, 20)};
fprintf ('results being saved in: %s_displacement\n \n', output_file)
sp_to_vtk (u, space, geometry, vtk_pts, sprintf ('%s_displacement', output_file), 'displacement')
%sp_to_vtk_stress (u, space, geometry, vtk_pts, problem_data.lambda_lame, ...
                  problem_data.mu_lame, sprintf ('%s_stress', output_file)); 

% Plot in Matlab
figure
def_geom = geo_deform (u, space, geometry);
nrbplot (def_geom.nurbs, [40 40 40], 'light', 'on')
axis equal tight
title ('Deformed configuration')