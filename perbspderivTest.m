%Test periodic basis spline derivative
U=[0 0 0 0.2 0.4 0.6 0.8 1 1 1];
u=0.08;
elem=5;
deriv=perbspderiv(U,u,elem,p);
evalPoint=periodicBasisEval(U,u,elem,p);
p=2;

hold on;

fplot(@(x)(periodicBasisEval(U, x, elem, p)), [0 0.999]);
fplot(@(x)(deriv*(x-u)+evalPoint), [0, 0.999]);