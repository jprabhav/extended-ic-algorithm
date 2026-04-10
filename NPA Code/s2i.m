function ind = s2i(D,S)
% S2I determines index by its decomposition in dimensions D
%
%   See also I2S.
%
ind = zeros(size(S,1),1);
for l=1:size(S,1)
    ind(l) = S(l,end);
    for k=(numel(D)-1):-1:1
        ind(l) = ind(l)+(S(l,k)-1)*prod(D(k+1:end));
    end
end
end