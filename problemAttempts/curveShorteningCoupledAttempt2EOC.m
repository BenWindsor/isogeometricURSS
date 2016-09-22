% Example problem from pozzi stinner

% Set up initial surface 
t=0;
xHandle=@(x)((1+0.25*exp(-t))*cos(2*pi*x));
yHandle=@(x)((1+0.25*exp(-t))*sin(2*pi*x));
degree=2;
elemNum=19;
N=[11 21 31 41 51];
T=1.0;
crvErrors=zeros(1,numel(N));
crvEOCs=zeros(1, numel(N)-1);
fieldErrors=zeros(1,numel(N));
fieldEOCs=zeros(1, numel(N)-1);
t=0.01;

for i=1:numel(N)
    
    % Set up time step etc.
    h=1/N(i);
    %delta=h*h;
    delta=h*h;
    sp=linspace(0, 1, N(i));
    
    %DEBUG
    fprintf('\n');
    fprintf(mat2str(sp));
    fprintf('\n');
    
    %Number of steps to get to time t
    steps=floor(t/delta);
    
    % Setup matrices to store the results
    storedxCoefs=zeros(N(i), steps);
    storedyCoefs=zeros(N(i), steps);
    storedFieldCoefs=zeros(N(i), steps);
    
    % Create initial geometry and test space
    prevCrv=periodicCurveInterpolate(N(i), degree, xHandle, yHandle);
    prevCrvCtrl=prevCrv.coefs;
    prevGeometry=geo_load(prevCrv);
    knots=prevGeometry.perbspline.knots;
    [qn, qw] = msh_set_quad_nodes(knots, msh_gauss_nodes(prevCrv.order));
    prevMsh=msh_cartesian(knots, qn, qw, prevGeometry);
    prevSpace=sp_perbsp(prevGeometry.perbspline, prevMsh);
    
    
    % Setup initial field values
    initialField=periodicCurveInterpolate(N(i), degree, @(x)(0.8));
    prevFieldCoefs=initialField.coefs';
    sourceTerm=@(x)(0);
    
    for j=1:steps
        % Calculate curve/surface stifness and mass matrices
        crvM = op_u_v_tp(prevSpace, prevSpace, prevMsh); %mu = 1
        crvA = op_gradu_gradv_tp(prevSpace, prevSpace, prevMsh); %epsilon = 1
        crvMat = (1/delta)*(crvM)+crvA;
        
        % Update surface
        prevxCoefs=prevCrvCtrl(1,:);
        prevyCoefs=prevCrvCtrl(2,:);
        
        %Create forcing terms for suface/curve update
        curve=perbspmak(prevCrvCtrl, knots);
        field=perbspmak(prevFieldCoefs', knots);
        forcing1=op_f_v_tp_param(prevSpace, prevMsh, @(x)(curveShorteningCoupledAttempt2Forcing1(curve, field, x)));
        forcing2=op_f_v_tp_param(prevSpace, prevMsh, @(x)(curveShorteningCoupledAttempt2Forcing2(curve, field, x)));
        
        crvRhs1=(1/delta)*crvM*prevxCoefs' + forcing1;
        crvRhs2=(1/delta)*crvM*prevyCoefs' + forcing2;
        
        newxCoefs=crvMat\crvRhs1;
        newyCoefs=crvMat\crvRhs2;
        storedxCoefs(:,steps)=newxCoefs;
        storedyCoefs(:,steps)=newyCoefs;
        
        newCrvCtrl=[newxCoefs'; newyCoefs'];
        
        % Create new crv and create new space on the new crv
        newCrv=perbspmak(newCrvCtrl, knots);
        newGeometry=geo_load(newCrv);
        [qn, qw] = msh_set_quad_nodes(knots, msh_gauss_nodes(newCrv.order));
        newMsh=msh_cartesian(knots, qn, qw, newGeometry);
        newSpace=sp_perbsp(newGeometry.perbspline, newMsh);
        
        % Update field
        fieldM=op_u_v_tp(newSpace, newSpace, newMsh);
        fieldA=op_gradu_gradv_tp(newSpace, newSpace, newMsh); %Q: on prev or new space??
        fieldMat = (1/delta)*fieldM + fieldA;
        
        %sourceTerm=@(x)(0); %Update source term in time
        
        fieldRHSM=op_u_v_tp(prevSpace, prevSpace, prevMsh);
        fieldRHS=(1/delta)*fieldRHSM*prevFieldCoefs;  %+fieldS;
        
        newFieldCoefs=fieldMat\fieldRHS;
        storedFieldCoefs(:,steps)=newFieldCoefs;
        
        % Update prev values to use new values
        prevCrvCtrl=newCrvCtrl;
        prevSpace=newSpace;
        prevMsh=newMsh;
        prevFieldCoefs=newFieldCoefs;
    end
    
    crvResults=[storedxCoefs(:,steps)'; storedyCoefs(:,steps)'];
    fieldResults=storedFieldCoefs(:,steps)';
    
    % Compare radii
    solutionHandle=@(x)(1+0.25*exp(-t));
    approxHandle=@(x)(norm(periodicSplineCurveEval(knots, x, degree, crvResults)));
    crvErrors(i)=maxDiff(solutionHandle, approxHandle, sp);
    %perbspplot(perbspmak(results, knots), 40);
    
    % Compare fields
    solutionHandle=@(x)(1.0/(1.0+0.25*exp(-t)));
    approxHandle=@(x)(periodicSplineCurveEval(knots, x, degree, fieldResults));
    fieldErrors(i)=maxDiff(solutionHandle, approxHandle, sp);
        
end

% Calculate the EOCs
for i=1:(numel(N)-1)
    numerator=log(crvErrors(i)/crvErrors(i+1));
    denominator=log((1/N(i))/(1/N(i+1)));
    crvEOCs(i)=numerator/denominator;
    
    numerator=log(fieldErrors(i)/fieldErrors(i+1));
    denominator=log((1/N(i))/(1/N(i+1)));
    fieldEOCs(i)=numerator/denominator;
end


% % Plot curves
% hold on;
% sp=linspace(0,1,100);
% for i=1:stepNum
%     % Approx curves
%     coefs=[storedxCoefs(:,i)'; storedyCoefs(:,i)'];
%     crv=perbspmak(coefs, knots);
%     perbspplot(crv, 30);
%     
%     % Actual curves
%     t=(i)*delta;
%     xPoints=(1+0.25*exp(-t))*cos(2*pi*sp);
%     yPoints=(1+0.25*exp(-t))*sin(2*pi*sp);
%     plot(xPoints, yPoints);
% end
% 
% figure
% hold on;
% title('Approx fields');
% for i=1:stepNum
%     field=perbspmak(storedFieldCoefs(:,i)', knots);
%     perbspplot(field, 1);
% end
% 
% % Print actual fields
% figure
% hold on;
% title('Actual fields');
% for i=1:stepNum
%     t = i*delta;
%     fplot(@(x)(1.0/(1.0+0.25*exp(-t))), [0 1]);
% end