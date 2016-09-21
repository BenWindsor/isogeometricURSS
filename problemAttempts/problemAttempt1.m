%create the curve
crv=nrbcirc(1);
%elevate degree giving better results, higher degree = larger M matrix
%higher degree seems to give more knots anyway
%knot insertion too??
crv=nrbdegelev(crv,2);

%load geometry as a unit circle and convert to correct format
geometry = geo_load(crv);
%load the knots
knots = geometry.nurbs.knots;
%set quadrature points
[qn, qw] = msh_set_quad_nodes(knots, msh_gauss_nodes(geometry.nurbs.order));
%create the mesh
msh = msh_cartesian(knots, qn, qw, geometry);
%create the basis space
space = sp_nurbs(geometry.nurbs, msh);

%set lhs matrix
delta=0.005; %time step, too big?? too small??
D=1; %diffusion constant?
M = op_u_v_tp(space, space, msh); %mu = 1
A = D*op_gradu_gradv_tp(space, space, msh); %epsilon = 1
mat = (1/delta)*(M)+A;

%set rhs
c_0=ones(sqrt(numel(M)),1); %initialise C_0 first value all ones
S=op_f_v_tp(space, msh, @(x,y)(x)); %wont allow constant source?
rhs=(1/delta)*M*c_0+S;

%set up array to store results
stepNumber=5;
results=zeros(numel(rhs),stepNumber);

%iterate the process
for i=1:stepNumber
    c_i=mat\rhs
    results(:,i)=c_i
    
    %update rhs
    rhs=(1/delta)*M*c_i+S
end

%plot data
nrbplot(crv, 100); %plot with 100 subdivisions
hold on;

%change intensity of colour depending on value of C at that point
data=results(:,1);
for i=1:sqrt(numel(M))
    step=sqrt(numel(M));
    
    %plot points equidistand around curve test
    %scatter(x,y,size,colour,fill)
    scatter(cos(i*2*pi*(1/step)), sin(i*2*pi*(1/step)), 25, [0,0,data(i)*1], 'filled');
end

    
    

