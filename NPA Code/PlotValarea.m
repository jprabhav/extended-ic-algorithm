function [q1,q2] = PlotValarea(d,nn,lvl,p0,p1,p2)
% PLOTVALEREA plots a region in the correlation space spanned by a mixture
% of three NS-boxes, specified below under the SDP constraints
% in the d2dd scenario
d = 3;
n_xy = [d, 2];
p0 = MargBox(p0,n_xy,d);
p1 = MargBox(p1,n_xy,d);
p2 = MargBox(p2,n_xy,d);
% the mixture is p0(1-q1-q2)+p1*q1+p2*q2

% SPD
M = NpaMatrix([d 2],[d d],lvl);
[Mp,Mr] = SplitMatrix(M);
q2 = zeros(1,nn);
q1 = zeros(1,nn);
for k=1:nn
    q2(k) = (k-1)/(nn-1);
    q1(k) = optPlotValarea(d,Mp,Mr,p0,p1,p2,q2(k));
end

end

function e = optPlotValarea(d,Mp,Mr,p0,p1,p2,q2)

nr = size(Mr,3);
X = sdpvar(1,1+nr);
G = Mr(:,:,1)*X(2);
for k=2:nr
    G = G+Mr(:,:,k)*X(k+1);
end
G = G+Mp(:,:,1);
% marginal probabilities
for k=1:(2+d)*(d-1)
    G = G+Mp(:,:,k+1)*(p0(k)*(1-X(1)-q2)+p1(k)*X(1)+p2(k)*q2);
end
% other probabilities
for k=1:2*d*(d-1)^2
    ind = (2+d)*(d-1)+k;
    G = G+Mp(:,:,1+ind)*(p0(ind)*(1-X(1)-q2)+p1(ind)*X(1)+p2(ind)*q2);
end

Con = [G>=0, X(1)>=0];
Obj = X(1);

sol = optimize(Con,Obj,sdpsettings('solver','mosek','verbose',0));
if sol.problem == 0
    e = value(Obj);
else
    %sol.info
    %yalmiperror(sol.problem)
    e = 0;
end

end
