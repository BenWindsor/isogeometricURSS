elemNum=19;
U=zeros(elemNum,1);
for i=1:elemNum
    U(i)=i/elemNum;
end
U=[0 0 0 U' 1 1]; %Only two 1's as one comes from the i/elemNum with i=elemNum

ux=zeros(elemNum,1); %+1 for the zero eval point
%Eval the exact funciton at 0/N, 1/N, ... (N-1)/N
for i=0:(elemNum-1)
    ux(i+1)=cos(2*pi*(i/elemNum));
end

%Eval the exact funciton at 0/N, 1/N, ... (N-1)/N
uy=zeros(elemNum,1); %+1 for the zero eval point
for i=0:(elemNum-1)
    uy(i+1)=sin(2*pi*(i/elemNum));
end

xMatrix=zeros(elemNum, elemNum);
yMatrix=zeros(elemNum, elemNum);
for i=1:elemNum
    for j=1:elemNum
        xMatrix(i, j)=periodicBasisEval(U, (i-1)/elemNum, j, 2);
        yMatrix(i, j)=periodicBasisEval(U, (i-1)/elemNum, j, 2);
    end
end

xCoords=xMatrix\ux;
yCoords=yMatrix\uy;

ctrl=[xCoords'; yCoords'];
perbspplot(perbspmak(ctrl, U),100);
hold on;
scatter(ux, uy);

