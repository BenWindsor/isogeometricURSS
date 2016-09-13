% OUTPUT: calculate the EOCs for 1 time step of the curve shortening flow

% Set up fixed surface parameters
xHandle=@(x)(cos(2*pi*x));
yHandle=@(x)(sin(2*pi*x));
degree=2;
N=[11 21 31 41 51];
T=1.0;
errors=zeros(1,numel(N));
EOCs=zeros(1, numel(N)-1);

t=0.1;
hold on;
sp=linspace(0,1,50);
xActual=sqrt(1-2*t)*cos(2*pi*sp);
yActual=sqrt(1-2*t)*sin(2*pi*sp);
plot(xActual, yActual);

for i=1:numel(N)
    
    % Set up time step etc.
    h=1/N(i);
    delta=h*h;
    sp=linspace(0, 1, N(i));
    
    %Number of steps to get to time t
    steps=floor(t/delta);
    
    %DEBUG
    fprintf(mat2str(sp));
    fprintf('\n');
    %END DEBUG
    
    
    
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
    newxCoefs=zeros(N(i), steps);
    newyCoefs=zeros(N(i), steps);
    
    for j=1:steps
        rhs1=(1/delta)*M*prevxCoefs' + 0; %+0 for no forcing term
        rhs2=(1/delta)*M*prevyCoefs' + 0;
        
        newxCoefs(:,j)=mat\rhs1;
        newyCoefs(:,j)=mat\rhs2;
        
        prevxCoefs=newxCoefs(:,j)';
        prevyCoefs=newyCoefs(:,j)';
        
        % Update M and A on new space
        coefs=[newxCoefs(:,j)'; newyCoefs(:,j)'];
        geometry=geo_load(perbspmak(coefs, knots));
        
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
    end
      
    results=[newxCoefs(:,steps)'; newyCoefs(:,steps)'];
    
    % Compare radii
    solutionHandle=@(x)(sqrt(1-2*t));
    approxHandle=@(x)(norm(periodicSplineCurveEval(knots, x, degree, results)));
    errors(i)=maxDiff(solutionHandle, approxHandle, sp);
    perbspplot(perbspmak(results, knots), 40);
    
end

% Calculate the EOCs
for i=1:numel(EOCs)
    numerator=log(errors(i)/errors(i+1));
    %numerator=log(errors(i+1)/errors(i));
    denominator=log((1/N(i))/(1/N(i+1)));
    %denominator=log((1/N(i+1))/(1/N(i)));
    EOCs(i)=numerator/denominator;
end

