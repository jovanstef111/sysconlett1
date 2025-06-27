function X=null1(A)

[nr,nc]=size(A);

if nc>0
nmin=min(nr,nc);
[U,S,V]=svd(A);
j=1;
for i=1:nmin
if abs(S(i,i))<1e-8 % beshe -8
    break
end
j=i+1;
end

X=V(:,j:nc);

else
    X=zeros(nc ,0);
end