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
N2 = 21; 

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
crv1=periodicCurveInterpolate(N1, degree, xHandle, yHandle);
crv2=periodicCurveInterpolate(N2, degree, xHandle, yHandle);
geometry1=geo_load(crv1);
geometry2=geo_load(crv2);
knots1=geometry1.perbspline.knots;
knots2=geometry2.perbspline.knots;

% Create mesh
[qn, qw] = msh_set_quad_nodes(knots1, msh_gauss_nodes(crv1.order));
msh1=msh_cartesian(knots1, qn, qw, geometry1);
[qn, qw] = msh_set_quad_nodes(knots2, msh_gauss_nodes(crv2.order));
msh2=msh_cartesian(knots2, qn, qw, geometry2);

% Create discrete space of periodic basis splines
space1=sp_perbsp(geometry1.perbspline, msh1);
space2=sp_perbsp(geometry2.perbspline, msh2);

% Assemble LHS
D=1; 
M1 = op_u_v_tp(space1, space1, msh1); %mu = 1
M1=M1(1:N1, 1:N1);
M2 = op_u_v_tp(space2, space2, msh2); %mu = 1
M2=M2(1:N2, 1:N2);
A1 = D*op_gradu_gradv_tp(space1, space1, msh1); %epsilon = 1
A1=A1(1:N1, 1:N1);
A2 = D*op_gradu_gradv_tp(space2, space2, msh2); %epsilon = 1
A2=A2(1:N2, 1:N2);
mat1 = (1/delta1)*(M1)+A1;
mat2 = (1/delta2)*(M2)+A2;

% Initial setup
curveC_01=periodicCurveInterpolate(N1, 2, @(x)(sin(6*pi*x)));
curveC_02=periodicCurveInterpolate(N2, 2, @(x)(sin(6*pi*x)));
c_01=curveC_01.coefs';
c_02=curveC_02.coefs';
t=0;
f1=@(x,y)((1+16*(t+1)*delta1)*cos(8*pi*inverseCircle(x,y))+(8-9*(t+1)*delta1)*sin(6*pi*inverseCircle(x,y)));
f2=@(x,y)((1+16*(t+1)*delta2)*cos(8*pi*inverseCircle(x,y))+(8-9*(t+1)*delta2)*sin(6*pi*inverseCircle(x,y)));
S1=op_f_v_tp(space1, msh1, f1); 
S1=S1(1:N1);
S2=op_f_v_tp(space2, msh2, f2);
S2=S2(1:N2);
rhs1=(1/delta1)*M1*c_01+S1;
rhs2=(1/delta2)*M2*c_02+S2;
rhs1=rhs1(1:N1);
rhs2=rhs2(1:N2);


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
    S=op_f_v_tp(space1, msh1, f );
    
    % resize S
    S=S(1:N1);
    
    % Solve
    rhs=(1/delta1)*M1*prevResults+S;
  
    % resize rhs
    rhs=rhs(1:N1);
    
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
    S=op_f_v_tp(space2, msh2, f );
    
    % resize S
    S=S(1:N2);
    
    % Solve
    rhs=(1/delta2)*M2*prevResults+S;
  
    % resize rhs
    rhs=rhs(1:N2);
    
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
    approxHandle=@(x)(periodicSplineCurveEval(knots1, x, degree, results1(:,i)'));
    maxVals1(i)=maxDiff(solutionHandle, approxHandle, sp1);
end

for i=1:(N2*N2+1)
    solutionHandle=@(x)((delta2*i)*cos(8*pi*x)+(1-(delta2*i))*sin(6*pi*x));
    approxHandle=@(x)(periodicSplineCurveEval(knots2, x, degree, results2(:,i)'));
    maxVals2(i)=maxDiff(solutionHandle, approxHandle, sp2);
end

EOCNumerator1=max(maxVals1);
EOCNumerator2=max(maxVals2);

EOCNumerator=log(EOCNumerator1/EOCNumerator2);

EOC=EOCNumerator/EOCDenominator;
disp(EOC);
