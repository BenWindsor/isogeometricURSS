function newArray = shiftArray( array, shift )
% Shift each element of the given array back shift number of times so for
% example shiftArrray([1 2 3 4], 2)= [ 3 4 1 2 ]

newArray=array;

for i=1:shift
    
    tempStart=newArray(1);
    
    for j=1:(numel(newArray)-1)
        newArray(j)=newArray(j+1);
    end
    
    newArray(end)=tempStart;

end

