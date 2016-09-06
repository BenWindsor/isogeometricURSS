% SOLUTION: c(x,t)=t*cos(8Pi*x)+(1-t)sin(6Pi*x)
% SOURCE: s(x,t)=(1+16t)cos(8Pi*x)+(8-9t)sin(6Pi*x)
% SURFACE: unit circle
% 
% Set up fixed surface parameters
xHandle=@(x)(cos(2*pi*x));
yHandle=@(x)(sin(2*pi*x));
degree=2;
elemNum=15;

% Final time in interval [0, T]
T=1; 

% Number of evaluation points
N1 = 11; 
N2 = 16; 

% Step size
h1=1/N1;  
h2=1/N2; 

% Eval points
sp1 = linspace(0,T,N1); 
sp2 = linspace(0,T,N2);

% Time step of order (step size)^2
delta1 = h1^2; 
delta2 = h2^2;

% Denominator of the EOC
EOCDenominator=log(h1/h2);

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
D=1; 
M = op_u_v_tp(space, space, msh); %mu = 1
M=M(1:elemNum, 1:elemNum);
A = D*op_gradu_gradv_tp(space, space, msh); %epsilon = 1
A=A(1:elemNum, 1:elemNum);
mat1 = (1/delta1)*(M)+A;
mat2 = (1/delta2)*(M)+A;

% Initial setup
curveC_0=periodicCurveInterpolate(elemNum, 2, @(x)(sin(6*pi*x)));
c_0=curveC_0.coefs';
t=0;
f1=@(x,y)((1+16*(t+1)*delta1)*cos(8*pi*inverseCircle(x,y))+(8-9*(t+1)*delta1)*sin(6*pi*inverseCircle(x,y)));
f2=@(x,y)((1+16*(t+1)*delta2)*cos(8*pi*inverseCircle(x,y))+(8-9*(t+1)*delta2)*sin(6*pi*inverseCircle(x,y)));
S1=op_f_v_tp(space, msh, f1); 
S1=S1(1:elemNum);
S2=op_f_v_tp(space, msh, f2);
S2=S2(1:elemNum);
rhs1=(1/delta1)*M*c_0+S1;
rhs2=(1/delta2)*M*c_0+S2;
rhs1=rhs1(1:elemNum);
rhs2=rhs2(1:elemNum);


% calculate results
initialResults1=mat1\rhs1;
initialResults2=mat2\rhs2;


% set up array to store first results
results1=zeros(numel(rhs1),N1*N1);

prevResults=initialResults1;

for step=1:(N1*N1)
    % update rhs
    t=delta1*(step+1);
    f=@(x,y)((1+16*t)*cos(8*pi*inverseCircle(x,y))+(8-9*t)*sin(6*pi*inverseCircle(x,y)));
    S=op_f_v_tp(space, msh, f );
    
    % resize S
    S=S(1:elemNum);
    
    % Solve
    rhs=(1/delta1)*M*prevResults+S;
  
    % resize rhs
    rhs=rhs(1:elemNum);
    
    % Solve equation
    c_i=mat1\rhs;
    results1(:,step)=c_i;
    
    % Store current time step for next rhs update
    prevResults=c_i;   
    
end

% set up array to store second results
results2=zeros(numel(rhs2),N2*N2);

prevResults=initialResults2;

for step=1:(N2*N2)
    % update rhs
    t=delta2*(step+1);
    f=@(x,y)((1+16*t)*cos(8*pi*inverseCircle(x,y))+(8-9*t)*sin(6*pi*inverseCircle(x,y)));
    S=op_f_v_tp(space, msh, f );
    
    % resize S
    S=S(1:elemNum);
    
    % Solve
    rhs=(1/delta2)*M*prevResults+S;
  
    % resize rhs
    rhs=rhs(1:elemNum);
    
    % Solve equation
    c_i=mat2\rhs;
    results2(:,step)=c_i;
    
    % Store current time step for next rhs update
    prevResults=c_i; 
    
end

maxVals1=zeros(N1*N1, 1);
maxVals2=zeros(N2*N2, 1);

% Add initial results for delta*1 time step
results1=[initialResults1 results1];
results2=[initialResults2 results2];

for i=1:(N1*N1+1)
    solutionHandle=@(x)((delta1*i)*cos(8*pi*x)+(1-(delta1*i))*sin(6*pi*x));
    approxHandle=@(x)(periodicSplineCurveEval(knots, x, degree, results1(:,i)'));
    maxVals1(i)=maxDiff(solutionHandle, approxHandle, sp1);
end

for i=1:(N2*N2+1)
    solutionHandle=@(x)((delta2*i)*cos(8*pi*x)+(1-(delta2*i))*sin(6*pi*x));
    approxHandle=@(x)(periodicSplineCurveEval(knots, x, degree, results2(:,i)'));
    maxVals2(i)=maxDiff(solutionHandle, approxHandle, sp2);
end

EOCNumerator1=max(maxVals1);
EOCNumerator2=max(maxVals2);

EOCNumerator=log(EOCNumerator1/EOCNumerator2);

EOC=EOCNumerator/EOCDenominator;
disp(EOC);
