% Curve shortening flow of the unit circle with no forcing term SINGLE STEP
% Given a set of delta time steps from big to small it advances up to the
% first step each time using the different deltas.

% Set up fixed surface parameters
xHandle=@(x)(cos(2*pi*x));
yHandle=@(x)(sin(2*pi*x));
degree=2;
elemNum=63;

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
delta=[0.005, 0.001, 0.0005 ];
xresults=zeros(elemNum, numel(delta));
yresults=zeros(elemNum, numel(delta));

for j=1:numel(delta)
    
    M = op_u_v_tp(space, space, msh); %mu = 1
    M=M(1:elemNum, 1:elemNum);
    A = op_gradu_gradv_tp(space, space, msh); %epsilon = 1
    A=A(1:elemNum, 1:elemNum);
    mat = (1/delta(j))*(M)+A;
    
    
    prevxCoefs=crv.coefs(1,:)';
    prevyCoefs=crv.coefs(2,:)';
    
    steps=floor(delta(1)/delta(j));
    
    newxCoefs=zeros(elemNum, steps);
    newyCoefs=zeros(elemNum, steps);
    
    % Calculate solutions
    for i=1:steps
        rhs1=(1/delta(j))*M*prevxCoefs + 0; %+0 for no forcing term
        rhs2=(1/delta(j))*M*prevyCoefs + 0;
        
        newxCoefs(:,i)=mat\rhs1;
        newyCoefs(:,i)=mat\rhs2;
        
        prevxCoefs=newxCoefs(:,i);
        prevyCoefs=newyCoefs(:,i);
        
        % Create surface and load geometry
        crv=perbspmak([prevxCoefs'; prevyCoefs'], knots);
        geometry=geo_load(crv);
        
        % Create mesh
        [qn, qw] = msh_set_quad_nodes(knots, msh_gauss_nodes(crv.order));
        msh=msh_cartesian(knots, qn, qw, geometry);
        
        % Create discrete space of periodic basis splines
        space=sp_perbsp(geometry.perbspline, msh);
        
        % Update M and A on new curve ready for next iteration
        M = op_u_v_tp(space, space, msh); %mu = 1
        M=M(1:elemNum, 1:elemNum);
        A = op_gradu_gradv_tp(space, space, msh); %epsilon = 1
        A=A(1:elemNum, 1:elemNum);
        mat = (1/delta(j))*(M)+A;  
    end
    
    xresults(:,j)=prevxCoefs;
    yresults(:,j)=prevyCoefs;
    
end

% Plot solutions
hold on;
for i=1:numel(delta)
    newCoefs=[xresults(:,i)'; yresults(:,i)'];
    %scatter(newxCoefs(i,:), newyCoefs(i,:), 'x');
    newCrv=perbspmak(newCoefs, knots);
    %perbspplot(crv, 100);
    perbspplot(newCrv, 100);
end


sp=linspace(0,1,300);
xActual=sqrt(1-2*delta(1))*cos(2*pi*sp);
yActual=sqrt(1-2*delta(1))*sin(2*pi*sp);
plot(xActual, yActual);

