crv=periodicCurveInterpolate(23,2,@(x)(cos(2*pi*x)), @(y)(sin(2*pi*y)));
perbspplot(crv,100);
hold on;
crv2=nrbcirc(1);
nrbplot(crv2,100);

figure;

geometry=geo_load(crv);
knots=geometry.perbspline.knots;
[qn, qw] = msh_set_quad_nodes(knots, msh_gauss_nodes(crv.order));
msh=msh_cartesian(knots, qn, qw, geometry);

numElems=numel(knots)-2*(crv.order-1)-1;

space=sp_perbsp(crv, msh);

delta=0.005; %time step, too big?? too small??
D=2*pi*2*pi; %diffusion constant?
M = op_u_v_tp(space, space, msh); %mu = 1
A = D*op_gradu_gradv_tp(space, space, msh); %epsilon = 1
mat = (1/delta)*(M)+A;

%As only numElems x numElems non zero elements cut off the remaining two cols and rows
mat=mat(1:numElems, 1:numElems);

Omega=3; %frequency in source term from our predefined c(x) = 2sin(2pix)+sin(2pi*omega*x) c hat function
c_0=zeros(sqrt(numel(M)),1); %initialise C_0 first value all ones
%S=op_f_v_tp(space, msh, @(x,y)(((-2/delta)-2*(2*pi)*(2*pi))*sin(2*pi*inverseCircle(x,y)) + ((1/delta)+(2*pi*Omega)*(2*pi*Omega))*sin(2*pi*Omega*inverseCircle(x,y)) )); %wont allow constant source?
S=op_f_v_tp_param(space, msh, @(x)(((-2/delta)-2*(2*pi)*(2*pi))*sin(2*pi*x) + ((1/delta)+(2*pi*Omega)*(2*pi*Omega))*sin(2*pi*Omega*x) )); %wont allow constant source?

rhs=(1/delta)*M*c_0+S;

%As only numElems non-zero elems cut off remaining two elements.
rhs=rhs(1:numElems);

%calculate results
results=mat\rhs;

field=perbspmak(results',knots);
perbspplot(field,10);

%Plot actual result
actualField=@(x)(-2*sin(2*pi*x)+sin(2*pi*Omega*x));
fplot(actualField, [0 1]);
