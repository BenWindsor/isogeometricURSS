% First attempt at cell motility simulation

% Set up initial surface
radius=1;
xHandle=@(x)(radius*cos(2*pi*x));
yHandle=@(x)(radius*sin(2*pi*x));
degree=2;
elemNum=59;

% Set up equation parameters
actualArea=pi*radius*radius;
D1=1.0;
D3=7.0;
gamma=25000; %Q: is it really x10^4 or should it be x10^-4
r1=0.01;
r3=0.013; % Q: did Rasa mean 13x10^-3 or 1.3x10^-3?
s1=0.0001;
s3=0.2;
s=0.0001;
b1=0.1;
b3=0.005;
k1=1.4;
tolerance=0.001; % Tolerance for lambda approx

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
steps=60;
alpha=0.001;

% Set up fields
initiala1=periodicCurveInterpolate(elemNum, degree, @(x)(0.5));
%initiala1=periodicCurveInterpolate(elemNum, degree, @(x)(sin(2*pi*x)));
preva1Coefs=initiala1.coefs';
% TEST
% prevalCoefs=0.5*ones(elemNum,1);
preva1Coefs(5)=0.504;
preva1Coefs(6)=0.502;
preva1Coefs(7)=0.508;
preva1Coefs(8)=0.512;
preva1Coefs(9)=0.508;
preva1Coefs(10)=0.504;

a1Vals=perbspeval(initiala1, knots); 
preva2=mean(a1Vals);
%preva2=0.5;
initiala3=periodicCurveInterpolate(elemNum, degree, @(x)(1.3));
preva3Coefs=initiala3.coefs';

% Set up matrices to store the new control points/coefs
xCoefsStore=zeros(elemNum, steps);
yCoefsStore=zeros(elemNum, steps);
a1CoefsStore=zeros(elemNum, steps);
a2Store=zeros(steps);
a3CoefsStore=zeros(elemNum, steps);

for step=1:steps
    disp(step);
    % Update cell surface
    
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
    forcingTerm=sparse([forcingTerm1; forcingTerm2]);
    %Q: negatie or postivie terms in lambda Funcs?
    lambdaTerm1=op_f_v_tp_param(prevSpace, prevMsh, @(x)(cellMotilityAttempt1lambdaFunc1(prevCrv, x)));
    lambdaTerm2=op_f_v_tp_param(prevSpace, prevMsh, @(x)(cellMotilityAttempt1lambdaFunc2(prevCrv, x)));
    incompleteRHS=(Mblock + Bblock)*prevCrvCoefs + forcingTerm;
   
    %Q: whats a good inital interval?
    % Bisection method on f(lambda)=|approxArea(lambda)-actualArea| on [a, b]
    a=-10;
    b=10;
    recalculate=1;
    while recalculate==1
        %Solve with the lhs lambda
        aLambdaTerm1 = a*lambdaTerm1;
        aLlambdaTerm2 = a*lambdaTerm2;
        aLambdaTerm=sparse([aLambdaTerm1; aLlambdaTerm2]);
        aCrvRhs = incompleteRHS + aLambdaTerm;
        aNewCrvCoefs=crvMat\aCrvRhs;
        % Calculate approx. area with new curve
        approxArea=aLambdaTerm1'*aNewCrvCoefs(1:elemNum) + aLlambdaTerm2'*aNewCrvCoefs((elemNum+1):2*elemNum);
        aDiff=actualArea-approxArea;
        
        %Solve with the rhs lambda
        bLambdaTerm1 = b*lambdaTerm1;
        bLlambdaTerm2 = b*lambdaTerm2;
        bLambdaTerm=sparse([bLambdaTerm1; bLlambdaTerm2]);
        bCrvRhs = incompleteRHS + bLambdaTerm;
        bNewCrvCoefs=crvMat\bCrvRhs;
        % Calculate approx. area with new curve
        approxArea=bLambdaTerm1'*bNewCrvCoefs(1:elemNum) + bLlambdaTerm2'*bNewCrvCoefs((elemNum+1):2*elemNum);
        bDiff=actualArea-approxArea;

        % Solve for mid lambda
        mid=(a+b)/2;
        midLambdaTerm1 = mid*lambdaTerm1;
        midLlambdaTerm2 = mid*lambdaTerm2;
        midLambdaTerm=sparse([midLambdaTerm1; midLlambdaTerm2]);
        midCrvRhs = incompleteRHS + midLambdaTerm;
        midNewCrvCoefs=crvMat\midCrvRhs;
        % % Calculate approx. area with new curve
        approxArea=midLambdaTerm1'*midNewCrvCoefs(1:elemNum) + midLlambdaTerm2'*midNewCrvCoefs((elemNum+1):2*elemNum);
        midDiff=actualArea-approxArea;
        
        % Test if its within tolerance
        if abs(midDiff)<=tolerance
            recalculate=0;
            lambda=mid;
            newCrvCoefs=midNewCrvCoefs;
        else
            if sign(midDiff)==sign(aDiff)
                a=mid;
            else
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
    newSpace=sp_perbsp(newGeometry.perbspline, newMsh);
    
    % Update fields on surface
    newM = op_u_v_tp(newSpace, newSpace, newMsh);
    %prevM = op_u_v_tp(prevSpace, prevSpace, prevMsh);
    prevM=crvM; %Re-use calculated values
    A = op_gradu_gradv_tp(newSpace, newSpace, newMsh);
    a1Source=op_f_v_tp_param(prevSpace, prevMsh, @(x)(cellMotilityAttempt1a1Source(preva1, preva2, preva3, gamma, r1, s, s1, s3, b1, x)));
    a3Source=op_f_v_tp_param(prevSpace, prevMsh, @(x)(cellMotilityAttempt1a3Source(preva1, preva3, gamma, r3, b3, x )));
    % Q: use prevCrv or newCrv here?
    redistTerm = op_vel_dot_gradu_v_tp_param(newSpace, newSpace, newMsh, @(x)(cellMotilityAttempt1AdvectionField(prevCrv, newCrv, delta, x)))';
    
    % % Update a1
    a1mat = sparse((1/delta)*newM + redistTerm + D1*A);
    a1rhs=sparse(((1/delta)*prevM)*preva1Coefs + a1Source);
    a1CoefsStore(:,step)=a1mat\a1rhs;
    
    % % Update a2
    a1Vals=perbspeval(preva1, knots); %Q: is this old a1 or new a1 to use here?
    a2Store(step)=mean(a1Vals);
    
    % % Update a3
    a3mat = sparse((1/delta)*newM + redistTerm + D3*A);
    a3rhs=sparse(((1/delta)*prevM)*preva3Coefs + a3Source);
    a3CoefsStore(:,step)=a3mat\a3rhs;
    
    % Update all the surface ctrl points
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
for i=1:10:steps
    crv=perbspmak([xCoefsStore(:,i)'; yCoefsStore(:,i)'], knots);
    perbspplot(crv, 80);
end

% % Print a1 curves
figure;
hold on;
title('a1 values');
for i=1:10:steps
    crv=perbspmak(a1CoefsStore(:,i)', knots);
    perbspplot(crv, 1);
end

% % Print a3 curves
figure;
hold on;
title('a3 values');
for i=1:10:steps
    crv=perbspmak(a3CoefsStore(:,i)', knots);
    perbspplot(crv, 1);
end