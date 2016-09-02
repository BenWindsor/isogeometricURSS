% Degree 2 test
U=[0 0 0 0.2 0.4 0.6 0.8 1.0 1.0 1.0];
p=2;
i=3;
span=findSpan(U,0.1);
hold on;
title('Degree 2 splines');
for i=1:7
    fplot(@(x)(basisSplineEval(U, x, i, p)), [0, 0.999]);
end
figure

title('Degree 2 periodic splines');
hold on;
for i=1:5
    fplot(@(x)(periodicBasisEval(U, x, i, p)), [0 0.9999]);
end
figure;

% Degree 3 test
U=[0 0 0 0 0.2 0.4 0.6 0.8 1.0 1.0 1.0 1.0];
p=3;
i=3;
span=findSpan(U,0.1);
hold on;
title('Degree 3 splines');
for i=1:8
    fplot(@(x)(basisSplineEval(U, x, i, p)), [0, 0.999]);
end
figure

hold on;
title('Degree 3 periodic splines');
for i=1:5
    fplot(@(x)(periodicBasisEval(U, x, i, p)), [0 0.9999]);
end