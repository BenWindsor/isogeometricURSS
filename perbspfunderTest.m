%A comparison of my perbspfunder.m with basisfunder.m to show they work the
%same way and give the same structure

%Set up curve and take points
U=[0 0 0 0.2 0.4 0.6 0.8 1 1 1];
p=2;
nders=0;
lin=linspace(0, 0.4, 40);
derivs=perbspfunder(U,lin,p,nders);
points=derivs(:,1,1);

%plot first set
hold on;
scatter(lin, points);
scatter(lin, derivs(:,1,2));
scatter(lin, derivs(:,1,3));

%set up values for the other basisfunder.m file for normal b splines
figure
mknots = length (knots)-1;
mcp    = -p - 1 + mknots;
ndof   = mcp + 1;

s     = findspan(mcp, p, lin, U);
derivs = basisfunder (s, p, lin, U, nders);

%plot second set
hold on;
scatter(lin, derivs(:,1,1));
scatter(lin, derivs(:,1,2));
scatter(lin, derivs(:,1,3));