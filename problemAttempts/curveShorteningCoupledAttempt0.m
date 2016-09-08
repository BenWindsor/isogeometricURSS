% Set up initial surface 
t=0;
xHandle=@(x)((1+0.25*exp(-t))*cos(2*pi*x));
yHandle=@(x)((1+0.25*exp(-t))*sin(2*pi*x));
degree=2;
elemNum=35;

% Create initial geometry and test space 
prevCrv=periodicCurveInterpolate(elemNum, degree, xHandle, yHandle);
prevCrvCtrl=prevCrv.coefs;
prevGeometry=geo_load(prevCrv);
knots=prevGeometry.perbspline.knots;
[qn, qw] = msh_set_quad_nodes(knots, msh_gauss_nodes(prevCrv.order));
prevMsh=msh_cartesian(knots, qn, qw, prevGeometry);
prevSpace=sp_perbsp(prevGeometry.perbspline, prevMsh);

% Set time step
delta=0.05;
stepNum=5;

% Calculate curve/surface stifness and mass matrices
crvM = op_u_v_tp(prevSpace, prevSpace, prevMsh); %mu = 1
crvM=crvM(1:elemNum, 1:elemNum);
crvA = op_gradu_gradv_tp(prevSpace, prevSpace, prevMsh); %epsilon = 1
crvA=crvA(1:elemNum, 1:elemNum);
crvMat = (1/delta)*(crvM)+crvA;

% Setup initial field values
prevFieldCoefs=0.8*ones(elemNum, 1);
sourceTerm=@(x,y)(0*x+0*y);

% Setup matrices to store the results
storedxCoefs=zeros(elemNum, stepNum);
storedyCoefs=zeros(elemNum, stepNum);
storedFieldCoefs=zeros(elemNum, stepNum);

for step=1:stepNum
    
    % Update surface
    t=step*delta;
    xHandle=@(x)((1+0.25*exp(-t))*cos(2*pi*x));
    yHandle=@(x)((1+0.25*exp(-t))*sin(2*pi*x));
    newCrv=periodicCurveInterpolate(elemNum, degree, xHandle, yHandle);
    
    newCrvCtrl=newCrv.coefs;
    storedxCoefs(:,step)=newCrvCtrl(1,:)';
    storedyCoefs(:,step)=newCrvCtrl(2,:)';
    
    % Create new crv and create new space on the new crv
    newGeometry=geo_load(newCrv);
    [qn, qw] = msh_set_quad_nodes(knots, msh_gauss_nodes(newCrv.order));
    newMsh=msh_cartesian(knots, qn, qw, newGeometry);
    newSpace=sp_perbsp(newGeometry.perbspline, newMsh);

    % Update field
    fieldM=op_u_v_tp(newSpace, newSpace, newMsh);
    fieldM=fieldM(1:elemNum, 1:elemNum);
    fieldA=op_gradu_gradv_tp(newSpace, newSpace, newMsh);
    fieldA=fieldA(1:elemNum, 1:elemNum);
    fieldMat = (1/delta)*fieldM + fieldA;
    
    sourceTerm=@(x,y)(0*x+0*y); %Update source term in time
    
    fieldRHSM=op_u_v_tp(prevSpace, prevSpace, prevMsh);
    fieldRHSM=fieldRHSM(1:elemNum, 1:elemNum);
    fieldS=op_f_v_tp(newSpace, newMsh, sourceTerm);
    fieldS=fieldS(1:elemNum);
    fieldRHS=(1/delta)*fieldRHSM*prevFieldCoefs+fieldS;
    
    newFieldCoefs=fieldMat\fieldRHS;
    storedFieldCoefs(:,step)=newFieldCoefs;
    
    % Update prev values to use new values
    prevSpace=newSpace;
    prevMsh=newMsh;
    prevFieldCoefs=newFieldCoefs;
    
end

% Print curves
hold on;
title('Curve evolution');
for i=1:stepNum
    ctrl=[storedxCoefs(:,i)'; storedyCoefs(:,i)'];
    perbspplot(perbspmak(ctrl, knots),30);
end

% Print approximated fields
figure
hold on;
title('Approximated fields');
for i=1:stepNum
    perbspplot(perbspmak(storedFieldCoefs(:,i)',knots),1);
end

% Print actual fields
figure
hold on;
title('Actual fields');
% Solve the ODE for actual solution
y0=0.8;
for i=1:stepNum
    tspan(i)=i*delta;
end
[t, y] = ode45(@(t,y) -y*(-0.25*exp(-t)), tspan, y0);
for i=1:stepNum
    fplot(@(x)(y(i)), [0 1]);
end