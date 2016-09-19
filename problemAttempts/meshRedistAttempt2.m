% A curve shortening flow but with the mesh point redistribution system,
% second approach in the rasa_project pdf
%Curve shortening flow of the unit circle with no forcing term SINGLE STEP

% Set up fixed surface parameters
% xHandle=@(x)(cos(2*pi*x));
% yHandle=@(x)(sin(2*pi*x));
xHandle=@(x)(nonUniformCircleXParam(x));
yHandle=@(x)(nonUniformCircleYParam(x));
degree=2;
elemNum=23;

% Create surface and load geometry
crv=periodicCurveInterpolate(elemNum, degree, xHandle, yHandle);
perbspplot(crv,100);
geometry=geo_load(crv);
knots=geometry.perbspline.knots;

% Create mesh
[qn, qw] = msh_set_quad_nodes(knots, msh_gauss_nodes(crv.order));
msh=msh_cartesian(knots, qn, qw, geometry);

% Create discrete space of periodic basis splines
space=sp_perbsp(geometry.perbspline, msh);

% Set time step
delta=0.01;
steps=50;

% Alpha param
alpha=0.01;
%sigma=0.1;

% Initial matrices
M = op_u_v_tp(space, space, msh); %mu = 1
M=M(1:elemNum, 1:elemNum);
A = op_gradu_gradv_tp(space, space, msh); %epsilon = 1
A=A(1:elemNum, 1:elemNum);
B11 = op_u_v_tp_param(space, space, msh, @(x)(meshRedistAttempt2func11B(crv, x)));
B11 = B11(1:elemNum, 1:elemNum);
B12 = op_u_v_tp_param(space, space, msh, @(x)(meshRedistAttempt2func12B(crv, x)));
B12 = B12(1:elemNum, 1:elemNum);
B21 = B12; %op_f_u_v_tp_param(space, msh, @(x)(meshRedistAttempt1func1(crv, x)));
B21 = B21(1:elemNum, 1:elemNum);
B22 = op_u_v_tp_param(space, space, msh, @(x)(meshRedistAttempt2func22B(crv, x)));
B22 = B22(1:elemNum, 1:elemNum);
Mblock = (alpha/delta) * [M sparse(elemNum,elemNum); sparse(elemNum,elemNum) M];
Bblock = ((1-alpha)/delta) * [B11 B12; B21 B22];
Ablock = [A sparse(elemNum,elemNum); sparse(elemNum,elemNum) A];
mat = sparse(Mblock + Bblock + Ablock);

% Put all coefs into one column 
prevxCoefs=crv.coefs(1,:);
prevyCoefs=crv.coefs(2,:);
prevCoefs=[prevxCoefs'; prevyCoefs'];


%newxCoefs=zeros(steps, elemNum);
%newyCoefs=zeros(steps, elemNum);

newCoefs=zeros(2*elemNum, steps);
% Calculate solutions
for i=1:steps
    rhs=(Mblock + Bblock)*prevCoefs;

    newCoefs(:,i)=mat\rhs;
    
    prevCoefs=newCoefs(:,i);
    
    % Create surface and load geometry
    crv=perbspmak([newCoefs(1:elemNum, i)'; newCoefs((elemNum+1):(2*elemNum),i)'], knots);
    geometry=geo_load(crv);
    
    % Create mesh
    [qn, qw] = msh_set_quad_nodes(knots, msh_gauss_nodes(crv.order));
    msh=msh_cartesian(knots, qn, qw, geometry);
    
    % Create discrete space of periodic basis splines
    space=sp_perbsp(geometry.perbspline, msh);
    
    % Update M, B and A on new curve ready for next iteration
    M = op_u_v_tp(space, space, msh); %mu = 1
    M=M(1:elemNum, 1:elemNum);
    A = op_gradu_gradv_tp(space, space, msh); %epsilon = 1
    A=A(1:elemNum, 1:elemNum);
    B11 = op_u_v_tp_param(space, space, msh, @(x)(meshRedistAttempt2func11B(crv, x)));
    B11 = B11(1:elemNum, 1:elemNum);
    B12 = op_u_v_tp_param(space, space, msh, @(x)(meshRedistAttempt2func12B(crv, x)));
    B12 = B12(1:elemNum, 1:elemNum);
    B21 = B12; %op_f_u_v_tp_param(space, msh, @(x)(meshRedistAttempt1func1(crv, x)));
    B21 = B21(1:elemNum, 1:elemNum);
    B22 = op_u_v_tp_param(space, space, msh, @(x)(meshRedistAttempt2func22B(crv, x)));
    B22 = B22(1:elemNum, 1:elemNum);
    Mblock = (alpha/delta) * [M sparse(elemNum,elemNum); sparse(elemNum,elemNum) M];
    Bblock = ((1-alpha)/delta) * [B11 B12; B21 B22];
    Ablock = [A sparse(elemNum,elemNum); sparse(elemNum,elemNum) A];
    mat = sparse(Mblock + Bblock + Ablock);
    
end

% Plot solutions
figure
hold on;
for i=1:1:steps
    currentCoefs=[newCoefs(1:elemNum,i)'; newCoefs((elemNum+1):(2*elemNum),i)'];
    scatter(newCoefs(1:elemNum,i), newCoefs((elemNum+1):(2*elemNum),i), 'x');
    newCrv=perbspmak(currentCoefs, knots);
    perbspplot(newCrv, 100);
end

% Plot actual
hold on;
for i=1:1:steps
    sp=linspace(0,1,50);
    xActual=sqrt(1-2*i*delta)*cos(2*pi*sp);
	yActual=sqrt(1-2*i*delta)*sin(2*pi*sp);
    plot(xActual, yActual);
end
