U=[0 0 0 0.2 0.4 0.6 0.8 1.0 1.0 1.0];
p=2;
i=3;

hold on;
for i=1:7
    fplot(@(x)(basisSplineEval(U, x, i, p)), [0, 0.999]);
end
% figure
% 
% hold on;
% for i=1:1
%     fplot(@(x)(periodicBasisEval(U, x, i, p)), [0 0.999]);
% end