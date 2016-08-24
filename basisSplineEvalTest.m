% Test the splineEval.m file
%example from pg55 of 'The NURBS Book' 
% U=[0,0,0,1,2,3,4,4,5,5,5];
% p=2;
% hold on;
% for i=1:8
%     fplot(@(x)(basisSplineEval(U,x,i,p)), [0,4.9999]);
% end

%example from pg67 of 'The NURBS Book'
U=[0,0,0,0,1,5,6,8,8,8,8];
p=3;
hold on;
for i=1:7
    fplot(@(x)(basisSplineEval(U,x,i,p)), [0,7.9999]);
end




