%Curve shortening flow of the unit circle with no forcing term SINGLE STEP

% Set up fixed surface parameters
xHandle=@(x)(cos(2*pi*x));
yHandle=@(x)(sin(2*pi*x));
degree=2;
elemNum=15;

% Create surface and load geometry
crv=periodicCurveInterpolate(elemNum, degree, xHandle, yHandle);
geometry=geo_load(crv);
knots=geometry.perbspline.knots;

% Create mesh
[qn, qw] = msh_set_quad_nodes(knots, msh_gauss_nodes(crv.order));
msh=msh_cartesian(knots, qn, qw, geometry);

% Create discrete space of periodic basis splines
space=sp_perbsp(geometry.perbspline, msh);

% Set time step
delta=0.05;

M = op_u_v_tp(space, space, msh); %mu = 1
M=M(1:elemNum, 1:elemNum);
A = op_gradu_gradv_tp(space, space, msh); %epsilon = 1
A=A(1:elemNum, 1:elemNum);
mat = (1/delta)*(M)+A;


prevxCoefs=crv.coefs(1,:);
prevyCoefs=crv.coefs(2,:);

steps=5;

newxCoefs=zeros(steps, elemNum);
newyCoefs=zeros(steps, elemNum);

% Calculate solutions
for i=1:steps
    rhs1=(1/delta)*M*prevxCoefs' + 0; %+0 for forcing term 
    rhs2=(1/delta)*M*prevyCoefs' + 0;

    newxCoefs(i,:)=mat\rhs1;
    newyCoefs(i,:)=mat\rhs2;
    
    prevxCoefs=newxCoefs(i,:);
    prevyCoefs=newyCoefs(i,:);
end

% Plot solutions
hold on;
for i=1:steps
    newCoefs=[newxCoefs(i,:); newyCoefs(i,:)];
    newCrv=perbspmak(newCoefs, knots);
    perbspplot(crv, 100);
    perbspplot(newCrv, 100);
end
