function result = op_f_u_v_tp_param( space, mesh, f)
% OP_F_V_TP: assemble the right-hand side vector r = [r(i)], with  r(i) = (f*u_i, v_i)
% taking values in the parametric domain.
% Derived from op_f_v_tp.m assemblng it column by column
% N.B. this is assuming the same space for u_i and v_i basis funcitons
% INPUT:
% space=space of basis functions 
% mesh=mesh of geometry
% f=function in inner product

% with coeff=identity function, this is just op_u_v_tp.m

% Get element number and set up knots
sp_univ=space.sp_univ;
spaceSize=size(sp_univ.connectivity);
elemNum=spaceSize(2);
U=sp_univ.knots;
p=sp_univ.degree;

result=zeros(elemNum, elemNum);
for i=1:elemNum
    % Create new function by combining f and the basis function for that
    % column
    newf=@(x)(f(x).*periodicBasisEval(U, x, i, p));
    temp=op_f_v_tp_param(space, mesh, newf);
    temp=temp(1:elemNum);
    result(:,i)=temp;
end

end

% Can test by checking:
% op_f_u_v_tp_param(space, mesh, @(x)(1+0*x))
% is the same as:
% op_u_v_tp(space, space, mesh)
