U=[0 0 0 1/7 2/7 3/7 4/7 5/7 6/7 1 1 1 ];
ctrl=[1.0000    0.2470   -0.6920   -1.1099   -0.6920    0.2470    1.0000; 0.4816    1.0821    0.8678    0.0000   -0.8678   -1.0821   -0.4816];
crv=perbspmak(ctrl, U);
perbspplot(crv,100);
hold on;
crv2=nrbcirc(1);
nrbplot(crv2,100);

numElems=numel(U)-2*(crv.order-1)-1;

figure;

geometry=geo_load(crv);
knots=geometry.perbspline.knots;
[qn qw] = msh_set_quad_nodes(knots, msh_gauss_nodes(crv.order));
msh=msh_cartesian(knots, qn, qw, geometry);

space=sp_perbsp(crv, msh);

delta=1.0; %time step, too big?? too small??
D=0.0; %2*pi*2*pi; %diffusion constant?
M = op_u_v_tp(space, space, msh); %mu = 1
A = D*op_gradu_gradv_tp(space, space, msh); %epsilon = 1
mat = (1/delta)*(M)+A;

%As only numElems x numElems non zero elements cut off the remaining two cols and rows
mat=mat(1:numElems, 1:numElems);

Omega=3; %frequency in source term from our predefined c(x) = 2sin(2pix)+sin(2pi*omega*x) c hat function
c_0=zeros(sqrt(numel(M)),1); %initialise C_0 first value all ones
S=op_f_v_tp(space, msh, @(x,y)(periodicBasisEval(U,inverseCircle(x,y),3,2)) ); %wont allow constant source?
rhs=(1/delta)*M*c_0+S;

%As only numElems non-zero elems cut off remaining two elements.
rhs=rhs(1:numElems);

%calculate results
results=mat\rhs;

field=perbspmak(results',U);
perbspplot(field,10);

%Plot actual result
actualField=@(x)(-2*sin(2*pi*x)+sin(2*pi*Omega*x));
fplot(actualField, [0 1]);
