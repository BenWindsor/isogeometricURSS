% SOLUTION: c(x,t)=t*cos(8Pi*x)+(1-t)sin(6Pi*x)
% SOURCE: s(x,t)=(1+16t)cos(8Pi*x)+(8-9t)sin(6Pi*x)
% SURFACE: unit circle

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
% QUESTION: no diffusion constant D?? set it to 1??
% QUESTION: what is the value of c_0?? all zeros

% Assemble LHS
D=1; 
M = op_u_v_tp(space, space, msh); %mu = 1
A = D*op_gradu_gradv_tp(space, space, msh); %epsilon = 1
mat = (1/delta)*(M)+A;

% Resize mat
M=M(1:elemNum, 1:elemNum);
mat=mat(1:elemNum, 1:elemNum);

% Get c_0
curveC_0=periodicCurveInterpolate(elemNum, 2, @(x)(sin(6*pi*x)));
% Assemble RHS
%c_0=ones(sqrt(numel(M)),1); %initialise C_0 first value all ones
c_0=curveC_0.coefs';
t=0;
% Provide source function for solution field on reference domain \hat{c}=tcos(8pix)+(1-t)sin(6pix)
%f=@(x,y)((1+16*(t+1)*delta)*cos(8*pi*inverseCircle(x,y))+(8-9*(t+1)*delta)*sin(6*pi*inverseCircle(x,y))); %QUESETION: start at t=0 here?
%S=op_f_v_tp(space, msh, f ); 
f=@(x)((1+16*(t+1)*delta)*cos(8*pi*x)+(8-9*(t+1)*delta)*sin(6*pi*x)); %QUESETION: start at t=0 here?
S=op_f_v_tp_param(space, msh, f ); 
S=S(1:elemNum);
rhs=(1/delta)*M*c_0+S;

% resize rhs
rhs=rhs(1:elemNum);

% calculate results
initialResults=mat\rhs;


% set up array to store results
stepNumber=5;
results=zeros(numel(rhs),stepNumber);

prevResults=initialResults;
for step=1:stepNumber
    % update rhs
    t=delta*(step+1);
    %f=@(x,y)((1+16*t)*cos(8*pi*inverseCircle(x,y))+(8-9*t)*sin(6*pi*inverseCircle(x,y)));
    %S=op_f_v_tp(space, msh, f );
    f=@(x)((1+16*t)*cos(8*pi*x)+(8-9*t)*sin(6*pi*x));
    S=op_f_v_tp_param(space, msh, f ); 
    % resize S
    S=S(1:elemNum);
    
    % Solve
    rhs=(1/delta)*M*prevResults+S;
  
    % resize rhs
    rhs=rhs(1:elemNum);
    
    % Solve equation
    c_i=mat\rhs;
    results(:,step)=c_i;
    
    % Store current time step for next rhs update
    prevResults=c_i;
    
    
end

% Print results
hold on;
% Inital solution for t=0
perbspplot(perbspmak(initialResults', knots), 1);
title('Approximated soltuions');
% Time step results
for i=1:stepNumber
    perbspplot(perbspmak(results(:,i)', knots), 1);
end
figure;

% Print some sample solutions
hold on;
title('actual solutions');
%Inital solution for t=0;
%fplot(@(x)((delta*0)*cos(8*pi*x)+(1-(delta*0))*sin(6*pi*x)), [0 1])
for i=1:(stepNumber+1)
    fplot(@(x)((delta*i)*cos(8*pi*x)+(1-(delta*i))*sin(6*pi*x)), [0 1]);
end