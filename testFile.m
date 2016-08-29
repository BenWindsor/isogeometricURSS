%General file for testing functions out

U=[0 0 0 0.2 0.4 0.6 0.8 1 1 1];
% p=2;
% hold on;
% for i=1:7
%     fplot(@(x)(basisSplineEval(U,x,i,p)), [0, 0.999])
% end
crv=perbspmak([0 1 0 -1 0; 1 0 -1 0 1], U);
geometry=geo_load(crv);
knots=geometry.perbspline.knots;
[qn qw] = msh_set_quad_nodes(knots, msh_gauss_nodes(3));
msh=msh_cartesian(knots, qn, qw, geometry);

space=sp_perbsp(crv, msh);
mat=op_gradu_gradv_tp(space, space, msh);
%S=op_f_v_tp(space, msh, @(x,y)(((-2/delta)-2*(2*pi)*(2*pi))*sin(2*pi*inverseCircle(x,y)) + ((1/delta)+(2*pi*Omega)*(2*pi*Omega))*sin(2*pi*Omega*inverseCircle(x,y)) )); %wont allow constant source?
test=op_u_v_tp(space, space, msh);
