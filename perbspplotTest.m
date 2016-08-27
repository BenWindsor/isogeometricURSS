%Test plotting periodic basis
% 2D ctrl 
U=[0 0 0 1 2 3 3 3];
ctrl=[0 0.25 0.75 1 1;
       0 1 0.6 0.2 0];
   
crv=perbspmak(ctrl, U);
perbspplot(crv,100);

%1D ctrl
figure
U=[0 0 0 0.2 0.4 0.6 0.8 1 1 1];
ctrl=[1 2 3 4 -3];
p=2;
crv=perbspmak(ctrl, U);
perbspplot(crv,100);