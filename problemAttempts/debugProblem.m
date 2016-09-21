%DESCRIPTION: a file to solve the specific problem c we found 

%create the curve
crv=nrbcirc(1);
%elevate degree giving better results, higher degree = larger M matrix
%higher degree seems to give more knots anyway
%knot insertion too??
crv=nrbdegelev(crv, 3);

%for knot refine example see solve_laplace_iso.m

%load geometry as a unit circle and convert to correct format
geometry = geo_load(crv);
%load the knots
knots = geometry.nurbs.knots;
%set quadrature points
[qn, qw] = msh_set_quad_nodes(knots, msh_gauss_nodes(geometry.nurbs.order));
%create the mesh
msh = msh_cartesian(knots, qn, qw, geometry);
%create the basis space
space = sp_nurbs(geometry.nurbs, msh);

test=op_u_v_tp(space, space, msh);
% %set lhs matrix
% delta=0.005; %time step, too big?? too small??
% D=2*pi*2*pi; %diffusion constant?
% M = op_u_v_tp(space, space, msh); %mu = 1
% A = D*op_gradu_gradv_tp(space, space, msh); %epsilon = 1
% mat = (1/delta)*(M)+A;
% 
% %set rhs
% Omega=3; %frequency in source term from our predefined c(x) = 2sin(2pix)+sin(2pi*omega*x) c hat function
% c_0=zeros(sqrt(numel(M)),1); %initialise C_0 first value all ones
% S=op_f_v_tp(space, msh, @(x,y)(((-2/delta)-2*(2*pi)*(2*pi))*sin(2*pi*inverseCircle(x,y)) + ((1/delta)+(2*pi*Omega)*(2*pi*Omega))*sin(2*pi*Omega*inverseCircle(x,y)) )); %wont allow constant source?
% rhs=(1/delta)*M*c_0+S;
% 
% %calculate results
% results=mat\rhs;
% 
% %create the field from the calculated results
% %should I use geometry.nurbs.knots or space.knots? the former works the
% %latter doesnt...?
% field=nrbmak(transpose(results), knots); %results needs to be a row so transpose
% 
% %sample some points and plot
% evalPoints=linspace(0,1,70);
% [points weights]=nrbeval(field, evalPoints);
% scatter(evalPoints, points(1,:));
% hold on;
% 
% %plot expected result as comparison
% actualField=@(x)(-2*sin(2*pi*x)+sin(2*pi*Omega*x))
% for i=1:numel(evalPoints)
%     actualPoints(i)=actualField(evalPoints(i))
% end
% scatter(evalPoints, actualPoints, 'filled');
% 
% %calculate and plot error too
% %error=actualPoints-points(1,:);
% %scatter(evalPoints, error, 'x');
% 
% %code used in examples
% %sp_to_vtk(results, space, geometry, 20*ones(msh.ndim,1), 'solution.vts', 'u');
% 
