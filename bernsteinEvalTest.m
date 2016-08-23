hold on;
for i=1:5
    fplot(@(x)(bernsteinEval(i, 5, x)), [0,1]);
end