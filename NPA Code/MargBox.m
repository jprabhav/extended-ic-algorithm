function p = MargBox(p,n_xy,d)
v = i2s([n_xy(1) d n_xy(2) d], 1:prod(n_xy)*d^2);
va = i2s([n_xy(1) d],1:n_xy(1)*d);
vb = i2s([n_xy(2) d],1:n_xy(2)*d);
% generate marginal probabilities
pa = zeros(1,n_xy(1)*d);
for x=1:n_xy(1)
    for a=1:d
        pa((x-1)*d+a) = sum(p(v(:,1)==x & v(:,2)==a & v(:,3)==1));
    end
end
pb = zeros(1,n_xy(2)*d);
for y=1:n_xy(2)
    for b=1:d
        pb((y-1)*d+b) = sum(p(v(:,3)==y & v(:,1)==1 & v(:,4)==b));
    end
end
p = p(v(:,2)~=d & v(:,4)~=d);
pa = pa(va(:,2)~=d);
pb = pb(vb(:,2)~=d);
p = [pa,pb,p]; 
end