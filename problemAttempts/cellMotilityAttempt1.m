% First attempt at cell motility simulation

% Set up initial surface
radius=1;
xHandle=@(x)(radius*cos(2*pi*x));
yHandle=@(x)(radius*sin(2*pi*x));
degree=1;
elemNum=29;

% Set up equation parameters
actualArea=pi*radius*radius;
D1=1.0;
D3=7.0;
gamma=0.00025;
r1=0.01;
r3=0.013; % Q: did Rasa mean 13x10^-3 or 1.3x10^-3?
s1=0.0001;
s3=0.2;
s=0.0001;
b1=0.1;
b3=0.005;
k1=1.4;
tolerance=0.001; % Tolerance for lambda approx

initiala1=periodicCurveInterpolate(elemNum, 2, @(x)(0.5));
preva1Coefs=initiala1.coefs';
% TEST
% prevalCoefs=0.5*ones(elemNum,1);
% preva1Coefs(1)=0.53;
% preva1Coefs(2)=0.55;
preva1Coefs(3)=0.51;
preva1Coefs(4)=0.52;
preva1Coefs(5)=0.53;
preva1Coefs(6)=0.52;
preva1Coefs(7)=0.51;


a1Vals=perbspeval(initiala1, knots); %Q: is this old a1 or new a1 to use here?
preva2=mean(a1Vals);
%preva2=0.5;
initiala3=periodicCurveInterpolate(elemNum, 2, @(x)(1.3));
preva3Coefs=initiala3.coefs';

% Create surface and load geometry
prevCrv=periodicCurveInterpolate(elemNum, degree, xHandle, yHandle);
prevGeometry=geo_load(prevCrv);
knots=prevGeometry.perbspline.knots;

% Put all ctrl points into one column ready to be solved for
prevxCoefs=prevCrv.coefs(1,:);
prevyCoefs=prevCrv.coefs(2,:);
prevCrvCoefs=[prevxCoefs'; prevyCoefs'];

% Create mesh and space
[qn, qw] = msh_set_quad_nodes(knots, msh_gauss_nodes(prevCrv.order));
prevMsh=msh_cartesian(knots, qn, qw, prevGeometry);
prevSpace=sp_perbsp(prevGeometry.perbspline, prevMsh);

% Set parameters delta=time step, alpha=mesh redist. coefficient
delta=0.0005;
steps=30;
alpha=0.2;

% Set up matrices to store the new control points/coefs
xCoefsStore=zeros(elemNum, steps);
yCoefsStore=zeros(elemNum, steps);
a1CoefsStore=zeros(elemNum, steps);
a2Store=zeros(steps);
a3CoefsStore=zeros(elemNum, steps);

for step=1:steps
    disp(step);
    % Update cell surface
    a=-10; % Initial interval [a,b] to search for lambda in
    b=10;
    recalculate=1;
    preva1=perbspmak(preva1Coefs', knots);
    preva3=perbspmak(preva3Coefs',knots);
    crvM=op_u_v_tp(prevSpace, prevSpace, prevMsh);
    crvA=op_gradu_gradv_tp(prevSpace, prevSpace, prevMsh);
    B11 = op_u_v_tp_param(prevSpace, prevSpace, prevMsh, @(x)(cellMotilityAttempt1Func11(prevCrv, x)));
    B11 = B11(1:elemNum, 1:elemNum);
    B12 = op_u_v_tp_param(prevSpace, prevSpace, prevMsh, @(x)(cellMotilityAttempt1Func12(prevCrv, x)));
    B12 = B12(1:elemNum, 1:elemNum);
    B21 = B12;
    B21 = B21(1:elemNum, 1:elemNum);
    B22 = op_u_v_tp_param(prevSpace, prevSpace, prevMsh, @(x)(cellMotilityAttempt1Func22(prevCrv, x)));
    B22 = B22(1:elemNum, 1:elemNum);
    Mblock = (alpha/delta) * [crvM sparse(elemNum,elemNum); sparse(elemNum,elemNum) crvM];
    Bblock = ((1-alpha)/delta) * [B11 B12; B21 B22];
    Ablock = [crvA sparse(elemNum,elemNum); sparse(elemNum,elemNum) crvA];
    crvMat = sparse(Mblock + Bblock + Ablock);
    forcingTerm1 = op_f_v_tp_param(prevSpace, prevMsh, @(x)(cellMotilityAttempt1forcing1(prevCrv, preva1, k1, x)));
    forcingTerm2 = op_f_v_tp_param(prevSpace, prevMsh, @(x)(cellMotilityAttempt1forcing2(prevCrv, preva1, k1, x)));
    forcingTerm=[forcingTerm1; forcingTerm2];
    
    while recalculate==1
        lambdaTerm1=op_f_v_tp_param(prevSpace, prevMsh, @(x)(cellMotilityAttempt1lambdaFunc1(prevCrv, x)));
        lambdaTerm2=op_f_v_tp_param(prevSpace, prevMsh, @(x)(cellMotilityAttempt1lambdaFunc2(prevCrv, x)));
        % Solve with the lhs lambda
        aLambdaTerm1 = a*lambdaTerm1;
        aLlambdaTerm2 = a*lambdaTerm2;
        aLambdaTerm=[aLambdaTerm1; aLlambdaTerm2];
        aCrvRhs = (Mblock + Bblock)*prevCrvCoefs + forcingTerm + aLambdaTerm;
        aNewCrvCoefs=crvMat\aCrvRhs;
        % % Calculate approx. area with new curve
        approxArea=aLambdaTerm1'*aNewCrvCoefs(1:elemNum) + aLlambdaTerm2'*aNewCrvCoefs((elemNum+1):2*elemNum);
        aDiff=actualArea-approxArea;
        
        % Solve with the rhs lambda
        bLambdaTerm1 = b*lambdaTerm1;
        bLlambdaTerm2 = b*lambdaTerm2;
        bLambdaTerm=[bLambdaTerm1; bLlambdaTerm2];
        bCrvRhs = (Mblock + Bblock)*prevCrvCoefs + forcingTerm + bLambdaTerm;
        bNewCrvCoefs=crvMat\aCrvRhs;
        % % Calculate approx. area with new curve
        approxArea=bLambdaTerm1'*bNewCrvCoefs(1:elemNum) + bLlambdaTerm2'*bNewCrvCoefs((elemNum+1):2*elemNum);
        bDiff=actualArea-approxArea;
        
        % Solve for mid lambda
        mid=(a+b)/2;
        midLambdaTerm1 = mid*lambdaTerm1;
        midLlambdaTerm2 = mid*lambdaTerm2;
        midLambdaTerm=[midLambdaTerm1; midLlambdaTerm2];
        midCrvRhs = (Mblock + Bblock)*prevCrvCoefs + forcingTerm + midLambdaTerm;
        midNewCrvCoefs=crvMat\midCrvRhs;
        % % Calculate approx. area with new curve
        approxArea=midLambdaTerm1'*midNewCrvCoefs(1:elemNum) + midLlambdaTerm2'*midNewCrvCoefs((elemNum+1):2*elemNum);
        midDiff=actualArea-approxArea;
        
        disp(midDiff);
        % Test if its within tolerance
        if abs(midDiff)<=tolerance
            recalculate=0;
            lambda=mid;
            newCrvCoefs=midNewCrvCoefs;
        else
            % get correct sign depending on if function
            % f=actualArea-approxArea is positive or negative
            % CHECK: monotonic inc or dec with lambda???
            if midDiff<=0
                a=mid;
                b=b;
            else
                a=a;
                b=mid;
            end
        end             
    end
    
    % Create the new space
    xCoefsStore(:,step)=newCrvCoefs(1:elemNum); 
    yCoefsStore(:,step)=newCrvCoefs((elemNum+1):(2*elemNum));
    newCrv=perbspmak([xCoefsStore(:,step)'; yCoefsStore(:,step)'], knots);
    newGeometry=geo_load(newCrv);
    [qn, qw] = msh_set_quad_nodes(knots, msh_gauss_nodes(newCrv.order));
    newMsh=msh_cartesian(knots, qn, qw, newGeometry);
    newSpace=sp_perbsp(prevGeometry.perbspline, newMsh);
    
    % Update fields on surface
    newM = op_u_v_tp(newSpace, newSpace, newMsh);
    %prevM = op_u_v_tp(prevSpace, prevSpace, prevMsh);
    prevM=crvM; %Re-use calculated values
    A = op_gradu_gradv_tp(newSpace, newSpace, newMsh);
    a1Source=op_f_v_tp_param(prevSpace, prevMsh, @(x)(cellMotilityAttempt1a1Source(preva1, preva2, preva3, gamma, r1, s, s1, s3, b1, x)));
    a3Source=op_f_v_tp_param(prevSpace, prevMsh, @(x)(cellMotilityAttempt1a3Source(preva1, preva3, gamma, r3, b3, x )));
    % Q: use prevCrv or newCrv here?
    redistTerm = op_vel_dot_gradu_v_tp_param(newSpace, newSpace, newMsh, @(x)(cellMotilityAttempt1AdvectionField(prevCrv, newCrv, delta, x)));
    mat = (1/delta)*newM + redistTerm + D1*A;
    
    % % Update a1
    a1rhs=((1/delta)*prevM)*preva1Coefs + a1Source;
    a1CoefsStore(:,step)=mat\a1rhs;
    
    % % Update a2
    a1Vals=perbspeval(preva1, knots); %Q: is this old a1 or new a1 to use here?
    a2Store(step)=mean(a1Vals);
    
    % % Update a3
    a3rhs=((1/delta)*prevM)*preva3Coefs + a3Source;
    a3CoefsStore(:,step)=mat\a3rhs;
    
    % Update all the surface values and flags
    recalculate=1;
    prevCrvCoefs=newCrvCoefs;
    prevMsh=newMsh;
    prevSpace=newSpace;
    preva1Coefs=a1CoefsStore(:,step);
    preva2=a2Store(step);
    preva3Coefs=a3CoefsStore(:,step);
    prevCrv=perbspmak([xCoefsStore(:,step)'; yCoefsStore(:,step)'], knots);
end

% Display results
% % Print Cell shapes
hold on;
title('Cell membrane');
for i=1:1:steps
    crv=perbspmak([xCoefsStore(:,i)'; yCoefsStore(:,i)'], knots);
    perbspplot(crv, 80);
end

% % Print a1 curves
figure;
hold on;
title('a1 values');
for i=1:1:steps
    crv=perbspmak(a1CoefsStore(:,i)', knots);
    perbspplot(crv, 1);
end

% % Print a3 curves
figure;
hold on;
title('a3 values');
for i=1:1:steps
    crv=perbspmak(a3CoefsStore(:,i)', knots);
    perbspplot(crv, 1);
end