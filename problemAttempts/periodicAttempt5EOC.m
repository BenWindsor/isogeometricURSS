% SOLUTION: c(x) = 2sin(2pix)+sin(2pi*omega*x)
% OUTPUT: calculate EOCs for the elliptic problem on stationary curve, so
%         just a single time step of a time independant field
% SURFACE: unit circle

% Set up fixed surface parameters
xHandle=@(x)(cos(2*pi*x));
yHandle=@(x)(sin(2*pi*x));
degree=2;
N=[11 21 31 41 51 61];
T=1;
%results=zeros(max(N),1);
errors=zeros(1,numel(N));
EOCs=zeros(1, numel(N)-1);

for i=1:numel(N)
    
    % Set up time step etc.
    h=1/N(i);
    sp=linspace(0, 1, N(i));
    delta=h*h;
    
    % Create surface and load geometry
    crv=periodicCurveInterpolate(N(i), degree, xHandle, yHandle);
    geometry=geo_load(crv);
    knots=geometry.perbspline.knots;
    
    % Create mesh
    [qn, qw] = msh_set_quad_nodes(knots, msh_gauss_nodes(crv.order));
    msh=msh_cartesian(knots, qn, qw, geometry);
    
    % Create discrete space of periodic basis splines
    space=sp_perbsp(geometry.perbspline, msh);
    
    % Assemble LHS
    D=2*pi*2*pi;
    M = op_u_v_tp(space, space, msh); %mu = 1
    A = D*op_gradu_gradv_tp(space, space, msh); %epsilon = 1
    mat = (1/delta)*(M)+A;
    % Resize mat
    M=M(1:N(i), 1:N(i));
    mat=mat(1:N(i), 1:N(i));
    
    Omega=3;
    c_0=zeros(sqrt(numel(M)),1);
    S=op_f_v_tp(space, msh, @(x,y)(((-2/delta)-2*(2*pi)*(2*pi))*sin(2*pi*inverseCircle(x,y)) + ((1/delta)+(2*pi*Omega)*(2*pi*Omega))*sin(2*pi*Omega*inverseCircle(x,y)) )); %wont allow constant source?
    S=S(1:N(i));
    rhs=(1/delta)*M*c_0+S;
    
    rhs=rhs(1:N(i));
    %calculate results
    results=mat\rhs;
    
    solutionHandle=@(x)(-2*sin(2*pi*x)+sin(2*pi*Omega*x));
    approxHandle=@(x)(periodicSplineCurveEval(knots, x, degree, results'));
    errors(i)=maxDiff(solutionHandle, approxHandle, sp);
    
    
%     field=perbspmak(transpose(results), knots); %results needs to be a row so transpose
%     sample some points and plot
%     perbspplot(field,1);
%     hold on;
%     plot expected result as comparison
%     fplot(@(x)(-2*sin(2*pi*x)+sin(2*pi*Omega*x)), [0 1])
%     figure
end

% Calculate the EOCs
for i=1:numel(EOCs)
    numerator=log(errors(i)/errors(i+1));
    denominator=log((1/N(i))/(1/N(i+1)));
    EOCs(i)=numerator/denominator;
end

