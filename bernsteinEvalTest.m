hold on;
% for i=1:5
%     fplot(@(x)(bernsteinEval(i, 1, x)), [0,1]);
% end

fplot(@(x)(bernsteinEval(2,2,x)), [0 1]);