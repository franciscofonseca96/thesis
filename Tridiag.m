function [TM]=Tridiag(n,d,u,b)

TM=spdiags([b*ones(n,1) d*ones(n,1) u*ones(n,1)], -1:1, n,n);

end