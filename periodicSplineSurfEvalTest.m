U=[0 0 0 0.2 0.4 0.6 0.8 1 1 1];
p=2;
divisions=40;
ctrl=zeros(5,5);
% Fill control points
for i=1:5
    for j=1:5
        ctrl(i,j)=i+j;
    end
end

xPoints=linspace(0, 1, divisions);
yPoints=linspace(0, 1, divisions);
zPoints=zeros(divisions,divisions);

for i=1:divisions
    for j=1:divisions
        point=cell(1,1);
        point{1, 1}=[xPoints(i), yPoints(j)];
        result=periodicSplineSurfEval(U, point, p, ctrl);
        zPoints(i,j)=result{1};
    end
end

hold on;
for i=1:divisions
    for j=1:divisions
        scatter3(xPoints(i), yPoints(j), zPoints(i,j));
    end 
end
        