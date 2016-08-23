% Test the splineEval.m file

% points=8000;
% space = linspace(0, 1, points);
% U = [0 0.1 0.2 0.3 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];
% secondDegree=zeros(1, points);
% secondDegree2=zeros(1, points);
% for i=1:(points-1)
%     secondDegree(i)=basisSplineEval(U, space(i), 2, 2);
%     secondDegree2(i)=basisSplineEval(U, space(i), 3, 1);
% end
% 
% scatter(space, secondDegree, 10); %make them size 10 so smaller
% hold on;
% scatter(space, secondDegree2, 10);

% U=[0,0,0.33,0.67,1,1];
% p=1;
% hold on;
% for i=1:4
%     %If you include 1 in interval it wont desplay final spline properly
%     fplot(@(x)(basisSplineEval(U,x,i,p)), [0 0.99999]);
% end


U = [0,0,0,0,0.33,0.67,1,1,1,1];
p=3;
hold on;
for i=1:6
    fplot(@(x)(basisSplineEval(U,x,i,p)), [0 0.9999]);
end




