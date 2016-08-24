function [ inserted, new ] = bezierKnots( U, p )
% Takes a knot vector and returns the knot vector and new knots
% needed for bezier extraction
% This means we take the internal knots and insert enough to have each
% internal knot with multiplicity p. 

% Count the occurences of each knot
singleValues=unique(U); %each entry of U only once
counts=histc(U,singleValues);
%extract only the data for the internal knots
internalCounts=counts(2:numel(counts)-1);
internalSingleValues=singleValues(2:numel(singleValues)-1);

%initialise outputs
inserted=0;
new=0;

%NOTE: preallocate for speed, should optimise
for i=1:numel(internalSingleValues)
    %create the new knot vector
    for j=1:p
        %add new val to end of new
        new(numel(new)+1)=internalSingleValues(i);
    end
    
    %create added knots
    for j=1:(p-internalCounts(i))
        inserted(numel(inserted)+1)=internalSingleValues(i);
    end

end

%strip the initial zero values
inserted=inserted(2:numel(inserted));
new=new(2:numel(new));

%Add the external knots back on by cutting the sections out of U and
%concatenating them with the new knot vector
new = cat(2, U(1:counts(1)), new, U(numel(U)-counts(numel(counts)):numel(U)));