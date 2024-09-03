function y=nanmean(x,dim)

p=isnan(x);
x(p)=0;
y=sum(x(:))/sum(~p(:));
