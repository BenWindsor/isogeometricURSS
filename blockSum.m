function result = blockSum( A, B )
% Take two matrices and produce the new matrix:
% (A 0)
% (0 B)

size1=size(A);
size2=size(B);

result=zeros(size1(1)+size2(1), size1(2)+size2(2));
resultSize=size(result);
result(1:size1(1), 1:size1(2))=A;
result((resultSize(1)-size2(1)+1):resultSize(1), (resultSize(2)-size2(2)+1):resultSize(2))=B;

end

