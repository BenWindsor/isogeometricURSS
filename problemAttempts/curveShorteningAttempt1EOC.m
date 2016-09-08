% OUTPUT: calculate the EOCs for 1 time step of the curve shortening flow

% Set up fixed surface parameters
xHandle=@(x)(cos(2*pi*x));
yHandle=@(x)(sin(2*pi*x));
degree=2;
N=[11 21 31 41 51 61];
T=1;
errors=zeros(1,numel(N));
EOCs=zeros(1, numel(N)-1);

for i=1:numel(N)
    
    % Set up time step etc.
    h=1/N(i);
    sp=linspace(0, 1, N(i));
    delta=h*h;
    
    crv=periodicCurveInterpolate(N(i), degree, xHandle, yHandle);
    geometry=geo_load(crv);
    knots=geometry.perbspline.knots;
    
    % Create mesh
    [qn, qw] = msh_set_quad_nodes(knots, msh_gauss_nodes(crv.order));
    msh=msh_cartesian(knots, qn, qw, geometry);
    
    % Create discrete space of periodic basis splines
    space=sp_perbsp(geometry.perbspline, msh);
    
    % Calculate the matrices
    M = op_u_v_tp(space, space, msh); %mu = 1
    M=M(1:N(i), 1:N(i));
    A = op_gradu_gradv_tp(space, space, msh); %epsilon = 1
    A=A(1:N(i), 1:N(i));
    mat = (1/delta)*(M)+A;
    
    
    prevxCoefs=crv.coefs(1,:);
    prevyCoefs=crv.coefs(2,:);
    
    % Calculate solutions
    rhs1=(1/delta)*M*prevxCoefs' + 0; %+0 for forcing term
    rhs2=(1/delta)*M*prevyCoefs' + 0;
    
    newxCoefs=mat\rhs1;
    newyCoefs=mat\rhs2;
    
    results=[newxCoefs'; newyCoefs'];
    
    % Compare radii
    solutionHandle=@(x)(sqrt(1-2*delta));
    approxHandle=@(x)(norm(periodicSplineCurveEval(knots, x, degree, results)));
    errors(i)=maxDiff(solutionHandle, approxHandle, sp);
    
end

% Calculate the EOCs
for i=1:numel(EOCs)
    numerator=log(errors(i)/errors(i+1));
    denominator=log((1/N(i))/(1/N(i+1)));
    EOCs(i)=numerator/denominator;
end

