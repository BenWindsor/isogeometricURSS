%Test plotting periodic basis
U=[0 0 0 1 2 3 3 3];
ctrl=[0 0.25 0.75 1 1;
       0 1 0.6 0.2 0];
   
crv=perbspmak(ctrl, U);
perbspplot(crv,100);