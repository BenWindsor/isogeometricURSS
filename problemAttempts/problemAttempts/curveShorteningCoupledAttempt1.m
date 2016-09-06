% Set up initial surface 
xHandle=@(x)(cos(2*pi*x));
yHandle=@(x)(sin(2*pi*x));
degree=2;
elemNum=15;

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
prevFieldCoefs=zeros(elemNum, 1);
sourceTerm=@(x)(x);

% Setup matrices to store the results
storedxCoefs=zeros(elemNum, stepNum);
storedyCoefs=zeros(elemNum, stepNum);
storedFieldCoefs=zeros(elemNum, stepNum);

for step=1:stepNum
    
    % Update surface
    prevxCoefs=prevCrvCtrl(1,:);
    prevyCoefs=prevCrvCtrl(2,:);
    
    crvRhs1=(1/delta)*crvM*prevxCoefs' + 0; %+0 for forcing term 
    crvRhs2=(1/delta)*crvM*prevyCoefs' + 0;

    newxCoefs=crvMat\crvRhs1;
    newyCoefs=crvMat\crvRhs2;
    storedxCoefs(:,step)=newxCoefs;
    storedyCoefs(:,step)=newyCoefs;
    
    newCrvCtrl=[newxCoefs'; newyCoefs'];
    
    % Create new crv and create new space on the new crv
    newCrv=perbspmak(newCrvCtrl, knots);
    newGeometry=geo_load(prevCrv);
    [qn, qw] = msh_set_quad_nodes(knots, msh_gauss_nodes(newCrv.order));
    newMsh=msh_cartesian(knots, qn, qw, newGeometry);
    newSpace=sp_perbsp(newGeometry.perbspline, newMsh);

    % Update field
    fieldM=op_u_v_tp(newSpace, newSpace, newMsh);
    fieldM=fieldM(1:elemNum, 1:elemNum);
    fieldA=op_gradu_gradv_tp(newSpace, newSpace, newMsh);
    fieldA=fieldA(1:elemNum, 1:elemNum);
    fieldMat = (1/delta)*fieldM + fieldA;
    
    sourceTerm=@(x)(x); %Update source term in time
    
    fieldRHSM=op_u_v_tp(prevSpace, prevSpace, prevMsh);
    fieldRHSM=fieldRHSM(1:elemNum, 1:elemNum);
    fieldS=op_f_v_tp(newSpace, newMsh, sourceTerm);
    fieldRHS=(1/delta)*fieldRHSM*prevFieldCoefs+fieldS;
    
    newFieldCoefs=fieldMat\fieldRHS;
    storedFieldCoefs(:,step)=newFieldCoefs;
    
    % Update prev values to use new values
    prevCrvCtrl=newCrvCtrl;
    prevSpace=newSpace;
    prevMsh=newMsh;
    prevFieldCoefs=newFieldCoefs;
    
end
