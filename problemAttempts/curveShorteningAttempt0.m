% SOLUTION: c(x,t)=t*cos(8Pi*x)+(1-t)sin(6Pi*x)
% SOURCE: s(x,t)=(1+16t)cos(8Pi*x)+(8-9t)sin(6Pi*x)
% SURFACE: unit circle

% Set up fixed surface parameters
xHandle=@(x)(cos(2*pi*x));
yHandle=@(x)(sin(2*pi*x));
degree=2;
elemNum=19;

% Create surface and load geometry
crv=periodicCurveInterpolate(elemNum, degree, xHandle, yHandle);
geometry=geo_load(crv);
knots=geometry.perbspline.knots;

% Create mesh
[qn, qw] = msh_set_quad_nodes(knots, msh_gauss_nodes(crv.order));
msh=msh_cartesian(knots, qn, qw, geometry);

% Create discrete space of periodic basis splines
space=sp_perbsp(geometry.perbspline, msh);

% Assemble LHS
delta=0.05;
D=1; 
M = op_u_v_tp(space, space, msh); %mu = 1
A = D*op_gradu_gradv_tp(space, space, msh); %epsilon = 1
mat = (1/delta)*(M)+A;
% Resize mat
M=M(1:elemNum, 1:elemNum);
mat=mat(1:elemNum, 1:elemNum);

% Get c_0
curveC_0=periodicCurveInterpolate(elemNum, degree, @(x)(sin(6*pi*x)));
% Assemble RHS
%c_0=ones(sqrt(numel(M)),1); %initialise C_0 first value all ones
c_0=curveC_0.coefs';
t=0;
% Provide source function for solution field on reference domain \hat{c}=tcos(8pix)+(1-t)sin(6pix)
f=@(x,y)((1+16*(t+1)*delta)*cos(8*pi*inverseCircle(x,y))+(8-9*(t+1)*delta)*sin(6*pi*inverseCircle(x,y))); %QUESETION: start at t=0 here?
S=op_f_v_tp(space, msh, f ); 
S=S(1:elemNum);
rhs=(1/delta)*M*c_0+S;

rhs=rhs(1:elemNum);
%calculate results
results=mat\rhs;

%create the field from the calculated results
%should I use geometry.nurbs.knots or space.knots? the former works the
%latter doesnt...?
field=perbspmak(transpose(results), knots); %results needs to be a row so transpose

%sample some points and plot
perbspplot(field,1);
hold on;

%plot expected result as comparison
fplot(@(x)((delta*1)*cos(8*pi*x)+(1-(delta*1))*sin(6*pi*x)), [0 1])

%calculate and plot error too
%error=actualPoints-points(1,:);
%scatter(evalPoints, error, 'x');

%code used in examples
%sp_to_vtk(results, space, geometry, 20*ones(msh.ndim,1), 'solution.vts', 'u');

