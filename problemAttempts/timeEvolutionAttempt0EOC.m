% SOLUTION: c(x,t)=t*cos(8Pi*x)+(1-t)sin(6Pi*x)
% SOURCE: s(x,t)=(1+16t)cos(8Pi*x)+(8-9t)sin(6Pi*x)
% SURFACE: unit circle

% Set up fixed surface parameters
xHandle=@(x)(cos(2*pi*x));
yHandle=@(x)(sin(2*pi*x));
degree=3;
N=[11 21 31];
T=1;
results=zeros(max(N),1);

for i=1:numel(N)
    
    % Set up time step etc.
    h=1/N(i);
    sp=linsapce(0, T, N(i));
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
    D=1;
    M = op_u_v_tp(space, space, msh); %mu = 1
    A = D*op_gradu_gradv_tp(space, space, msh); %epsilon = 1
    mat = (1/delta)*(M)+A;
    % Resize mat
    M=M(1:N(i), 1:N(i));
    mat=mat(1:N(i), 1:N(i));
    
    % Get c_0
    curveC_0=periodicCurveInterpolate(elemNum, degree, @(x)(sin(6*pi*x)));
    % Assemble RHS
    %c_0=ones(sqrt(numel(M)),1); %initialise C_0 first value all ones
    c_0=curveC_0.coefs';
    t=0;
    % Provide source function for solution field on reference domain \hat{c}=tcos(8pix)+(1-t)sin(6pix)
    f=@(x,y)((1+16*(t+1)*delta)*cos(8*pi*inverseCircle(x,y))+(8-9*(t+1)*delta)*sin(6*pi*inverseCircle(x,y))); %QUESETION: start at t=0 here?
    S=op_f_v_tp(space, msh, f );
    S=S(1:N(i));
    rhs=(1/delta)*M*c_0+S;
    
    rhs=rhs(1:N(i));
    %calculate results
    results(i)=mat\rhs;

end

field=perbspmak(transpose(results), knots); %results needs to be a row so transpose

%sample some points and plot
perbspplot(field,1);
hold on;

%plot expected result as comparison
fplot(@(x)((delta*1)*cos(8*pi*x)+(1-(delta*1))*sin(6*pi*x)), [0 1])

