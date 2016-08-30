U=[0         0         0    0.0526    0.1053    0.1579    0.2105    0.2632    0.3158    0.3684    0.4211    0.4737    0.5263 0.5789    0.6316    0.6842    0.7368    0.7895    0.8421    0.8947    0.9474    1.0000    1.0000    1.0000];
ctrl=[0.8916    0.6866    0.4072    0.0837   -0.2489   -0.5545   -0.8001   -0.9589   -1.0138   -0.9589   -0.8001   -0.5545   -0.2489   0.0837    0.4072    0.6866    0.8916    1.0000    1.0000;
    0.4825    0.7459    0.9284    1.0104    0.9828    0.8487    0.6227    0.3292    0.0000   -0.3292   -0.6227   -0.8487   -0.9828 -1.0104   -0.9284   -0.7459   -0.4825   -0.1669    0.1669];
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

delta=0.005; %time step, too big?? too small??
D=2*pi*2*pi; %diffusion constant?
M = op_u_v_tp(space, space, msh); %mu = 1
A = D*op_gradu_gradv_tp(space, space, msh); %epsilon = 1
mat = (1/delta)*(M)+A;

%As only numElems x numElems non zero elements cut off the remaining two cols and rows
mat=mat(1:numElems, 1:numElems);

Omega=3; %frequency in source term from our predefined c(x) = 2sin(2pix)+sin(2pi*omega*x) c hat function
c_0=zeros(sqrt(numel(M)),1); %initialise C_0 first value all ones
S=op_f_v_tp(space, msh, @(x,y)(((-2/delta)-2*(2*pi)*(2*pi))*sin(2*pi*inverseCircle(x,y)) + ((1/delta)+(2*pi*Omega)*(2*pi*Omega))*sin(2*pi*Omega*inverseCircle(x,y)) )); %wont allow constant source?
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
