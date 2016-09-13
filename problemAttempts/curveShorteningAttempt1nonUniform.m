%Curve shortening flow of the unit circle with no forcing term SINGLE STEP

% Set up fixed surface parameters
xHandle=@(x)(nonUniformCircleXParam(x));
yHandle=@(x)(nonUniformCircleYParam(x));
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
    rhs1=(1/delta)*M*prevxCoefs' + 0; %+0 for no forcing term 
    rhs2=(1/delta)*M*prevyCoefs' + 0;

    newxCoefs(i,:)=mat\rhs1;
    newyCoefs(i,:)=mat\rhs2;
    
    prevxCoefs=newxCoefs(i,:);
    prevyCoefs=newyCoefs(i,:);
    
    
    % Create surface and load geometry
    crv=perbspmak([prevxCoefs; prevyCoefs], knots);
    geometry=geo_load(crv);
    knots=geometry.perbspline.knots;
    
    % Create mesh
    [qn, qw] = msh_set_quad_nodes(knots, msh_gauss_nodes(crv.order));
    msh=msh_cartesian(knots, qn, qw, geometry);
    
    % Create discrete space of periodic basis splines
    space=sp_perbsp(geometry.perbspline, msh);
    
    % Update M and A on new curve ready for next iteration
    M = op_u_v_tp(space, space, msh); %mu = 1
    M=M(1:elemNum, 1:elemNum);
    A = op_gradu_gradv_tp(space, space, msh); %epsilon = 1
    A=A(1:elemNum, 1:elemNum);
    mat = (1/delta)*(M)+A;

end

% Plot solutions
hold on;
for i=1:steps
    newCoefs=[newxCoefs(i,:); newyCoefs(i,:)];
    scatter(newxCoefs(i,:), newyCoefs(i,:), 'x');
    newCrv=perbspmak(newCoefs, knots);
    %perbspplot(crv, 100);
    perbspplot(newCrv, 100);
end

% Plot actual
hold on;
for i=1:steps
    sp=linspace(0,1,50);
    xActual=sqrt(1-2*i*delta)*cos(2*pi*sp);
	yActual=sqrt(1-2*i*delta)*sin(2*pi*sp);
    plot(xActual, yActual);
end
