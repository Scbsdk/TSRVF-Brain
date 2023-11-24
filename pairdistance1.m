function D = pairdistance1(X,Y)
n1 = size(X,3);
n2 = size(Y,3);
D= zeros(n1,n2);

for i = 1:n1
    for j = 1:n2
        D(i,j)=DynamicQ(X(:,:,i),Y(:,:,j));
     end
end