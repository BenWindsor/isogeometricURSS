%Test periodic basis spline derivative
U=[0 0 0 0.2 0.4 0.6 0.8 1 1 1];
u=0.08;
elem=5;
deriv=perbspderiv(U,u,elem,p);
evalPoint=periodicBasisEval(U,u,elem,p);
p=2;

hold on;

%plot the basis function in question
fplot(@(x)(periodicBasisEval(U, x, elem, p)), [0 0.999]);
%plot the derivative tangent line to the point chosen above as u
fplot(@(x)(deriv*(x-u)+evalPoint), [0, 0.999]);