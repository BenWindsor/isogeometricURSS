function perbsp = periodicCurveInterpolate( elemNum, p, varargin )
% Interpolate curve given by function handles varargin with elemNum
% elements in the knot vector and degree p. Interpolation done with
% periodic basis splines.
% INPUTS
% elemNum=number of desired elements in knot span
% p=degree of periodic basis splines to use
% varargin=either a single funciton handle for a 1d curve or two function
% handles for a 2d parametric curve to approximate.

%N.B. Doesnt seem to do well with 1D curves, check
%periodicCurveInterpolateTest.m to see this in the commented out functions

%Account for the fact that even elements wont work past 20:
if (mod(elemNum,2)==0 && elemNum>19)
    error('please use an odd number of elements to avoid singular matrices');
end

%Deal with different inputs
if(numel(varargin)==1)
    xHandle=varargin{1};   
elseif(numel(varargin)==2)
    xHandle=varargin{1};
    yHandle=varargin{2};
else
    error('please give 1 or 2 function handles');
end

if(p~=2)
    error('Can only handle degree 2 at the moment');
end

U=zeros(elemNum,1);
startZeros=zeros(p+1,1);
endOnes=ones(p,1); %only p as one will come from the iteration below
for i=1:elemNum
    U(i)=i/elemNum;
end
U=[startZeros' U' endOnes']; %uniform open knot vector generated

%1D case
if(numel(varargin)==1)
    ux=zeros(elemNum,1); %+1 for the zero eval point
    %Eval the exact funciton at 0/N, 1/N, ... (N-1)/N
    for i=0:(elemNum-1)
        ux(i+1)=xHandle(i/elemNum);
    end
    
    xMatrix=zeros(elemNum, elemNum);
    for i=1:elemNum
        for j=1:elemNum
            xMatrix(i, j)=periodicBasisEval(U, (i-1)/elemNum, j, 2);
        end
    end
    
    xCoords=xMatrix\ux;
    
    ctrl=xCoords';
    
    perbsp=perbspmak(ctrl, U);
elseif (numel(varargin)==2)
     ux=zeros(elemNum,1); %+1 for the zero eval point
     uy=zeros(elemNum,1);
    %Eval the exact funciton at 0/N, 1/N, ... (N-1)/N
    for i=0:(elemNum-1)
        ux(i+1)=xHandle(i/elemNum);
        uy(i+1)=yHandle(i/elemNum);
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
    
    perbsp=perbspmak(ctrl, U);

end

