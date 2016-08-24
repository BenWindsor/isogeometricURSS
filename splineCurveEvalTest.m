% U=[0,0,0,0,1/4,1/2,3/4,1,1,1,1];
% p=3;
% ctrl=[0 1 1 2 3 2 1;
%       0 0 1 1 0 -2 -1];
% func=@(x)(splineCurveEval(U,x,p,ctrl));
% points=func(linspace(0,0.9999,200));
% scatter(points(1,:), points(2,:));

% U=[0,0,0,0,1,1,1,1];
% p=3;
% ctrl=[0 1 2 1;
%       0 2 3 0];
% func=@(x)(splineCurveEval(U,x,p,ctrl));
% points=func(linspace(0,0.999,200));
% scatter(points(1,:), points(2,:));
% hold on;
% for i=4:4
%     fplot(@(x)(basisSplineEval(U,x,i,p)),[0 0.999]);
%     
% end
%     

%example roughly from pg85 of 'The NURBS Book'
U=[0,0,0,1/5,2/5,3/5,4/5,1,1,1];
p=2;
ctrl=[0 -0.5 3 2 0.5 4 3.5;
      0 2 2.5 0 -2 -2 0];
func=@(x)(splineCurveEval(U,x,p,ctrl));
points=func(linspace(0,0.9999,200));
scatter(points(1,:), points(2,:));

