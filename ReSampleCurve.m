function Xn = ReSampleCurve(X,N)

    [n,T] = size(X);
    for r = 2:T
       del(r-1) = norm(X(:,r) - X(:,r-1));
    end
    cumdel = cumsum(del)/sum(del);
    cumdel=[0,cumdel];
    
    newdel = [1:N]/(N);
    for j=1:n
        Xn(j,:) = interp1(cumdel,X(j,1:T),newdel,'linear');
    end