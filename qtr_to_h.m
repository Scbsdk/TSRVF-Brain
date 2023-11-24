function h= qtr_to_h(q)
[m,n,T]=size(q);
for k=1:T
    for i=1:n-1
        gam=DynamicProgrammingQ(q(:,i,k)',q(:,i+1,k)',0,0);
        gam = invertGamma(gam);
        gam = (gam-gam(1))/(gam(end)-gam(1));
        temp1=Group_Action_by_Gamma_Coord(q(:,i+1,k)',gam).*sqrt(gradient(gam,1/(m-1)));
        alphadot(:,i,k)=temp1-q(:,i,k)';
        L(i,k)=sqrt(norm(alphadot(:,i,k)));
        if L(i,k)>0.0001
            h(:,i,k)=alphadot(:,i,k)/L(i,k);
        else
            h(:,i,k)=alphadot(:,i,k)*0.0001;
        end
    end
    h(:,n,k)=h(:,n-1,k);
end
