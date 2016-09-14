% Example problem from pozzi stinner

% Set up initial surface 
t=0;
xHandle=@(x)((1+0.25*exp(-t))*cos(2*pi*x));
yHandle=@(x)((1+0.25*exp(-t))*sin(2*pi*x));
degree=2;
elemNum=19;

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

% Setup initial field values
initialField=periodicCurveInterpolate(elemNum, 2, @(x)(0.8));
prevFieldCoefs=initialField.coefs';
sourceTerm=@(x)(0);

% Setup matrices to store the results
storedxCoefs=zeros(elemNum, stepNum);
storedyCoefs=zeros(elemNum, stepNum);
storedFieldCoefs=zeros(elemNum, stepNum);

for step=1:stepNum
    
    % Calculate curve/surface stifness and mass matrices
    crvM = op_u_v_tp(prevSpace, prevSpace, prevMsh); %mu = 1
    crvM=crvM(1:elemNum, 1:elemNum);
    crvA = op_gradu_gradv_tp(prevSpace, prevSpace, prevMsh); %epsilon = 1
    crvA=crvA(1:elemNum, 1:elemNum);
    crvMat = (1/delta)*(crvM)+crvA;
  
    % Update surface
    prevxCoefs=prevCrvCtrl(1,:);
    prevyCoefs=prevCrvCtrl(2,:);
    
    %Create forcing terms for suface/curve update
    curve=perbspmak(prevCrvCtrl, knots);
    field=perbspmak(prevFieldCoefs', knots);
    forcing1=op_f_v_tp_param(prevSpace, prevMsh, @(x)(curveShorteningCoupledAttempt2Forcing1(curve, field, x)));
    forcing1=forcing1(1:elemNum);
    forcing2=op_f_v_tp_param(prevSpace, prevMsh, @(x)(curveShorteningCoupledAttempt2Forcing2(curve, field, x)));
    forcing2=forcing2(1:elemNum);
    
    crvRhs1=(1/delta)*crvM*prevxCoefs' + forcing1; 
    crvRhs2=(1/delta)*crvM*prevyCoefs' + forcing2;

    newxCoefs=crvMat\crvRhs1;
    newyCoefs=crvMat\crvRhs2;
    storedxCoefs(:,step)=newxCoefs;
    storedyCoefs(:,step)=newyCoefs;
    
    newCrvCtrl=[newxCoefs'; newyCoefs'];
    
    % Create new crv and create new space on the new crv
    newCrv=perbspmak(newCrvCtrl, knots);
    newGeometry=geo_load(newCrv);
    [qn, qw] = msh_set_quad_nodes(knots, msh_gauss_nodes(newCrv.order));
    newMsh=msh_cartesian(knots, qn, qw, newGeometry);
    newSpace=sp_perbsp(newGeometry.perbspline, newMsh);
    
    % Update field
    fieldM=op_u_v_tp(newSpace, newSpace, newMsh);
    fieldM=fieldM(1:elemNum, 1:elemNum);
    fieldA=op_gradu_gradv_tp(newSpace, newSpace, newMsh); %Q: on prev or new space??
    fieldA=fieldA(1:elemNum, 1:elemNum);
    fieldMat = (1/delta)*fieldM + fieldA;
    
    %sourceTerm=@(x)(0); %Update source term in time
    
    fieldRHSM=op_u_v_tp(prevSpace, prevSpace, prevMsh);
    fieldRHSM=fieldRHSM(1:elemNum, 1:elemNum);
    %fieldS=op_f_v_tp(newSpace, newMsh, sourceTerm);
    %fieldS=fieldS(1:elemNum);
    fieldRHS=(1/delta)*fieldRHSM*prevFieldCoefs;  %+fieldS;
    
    newFieldCoefs=fieldMat\fieldRHS;
    storedFieldCoefs(:,step)=newFieldCoefs;
    
    % Update prev values to use new values
    prevCrvCtrl=newCrvCtrl;
    prevSpace=newSpace;
    prevMsh=newMsh;
    prevFieldCoefs=newFieldCoefs;
    
end

% Plot curves
hold on;
sp=linspace(0,1,50);
for i=1:stepNum
    % Approx curves
    coefs=[storedxCoefs(:,i)'; storedyCoefs(:,i)'];
    crv=perbspmak(coefs, knots);
    perbspplot(crv, 30);
    
    % Actual curves
    t=i*delta;
    xPoints=(1+0.25*exp(-t))*cos(2*pi*sp);
    yPoints=(1+0.25*exp(-t))*sin(2*pi*sp);
    plot(xPoints, yPoints);
end

figure
hold on;
for i=1:stepNum
    field=perbspmak(storedFieldCoefs(:,i)', knots);
    perbspplot(field, 1);
end

% Print actual fields
figure
hold on;
title('Actual fields');
for i=1:stepNum
    t = i*delta;
    fplot(@(x)(1.0/(1.0+0.25*exp(-t))), [0 1]);
end