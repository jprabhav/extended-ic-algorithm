function [M,T] = NpaMatrix(n_xy,n_ab,lvl)
% NPAMATRIX generate NPA matix for n_a effects of Alice and n_b of Bob
%   lvl -- level of NPA hierachy
%
% generate a 2-column cell array all the effects of Alice and Bob
% for the level lvl
%
n_a = n_xy(1)*(n_ab(1)-1); % number of effects of Alice
n_b = n_xy(2)*(n_ab(2)-1); % of Bob
T{1,1} = 0;
T{1,2} = 0;
for k=1:n_a
    T{k+1,1} = k;
    T{k+1,2} = 0;
end
for k=1:n_b
    T{n_a+1+k,1} = 0;
    T{n_a+1+k,2} = k;
end
T0 = T(2:end,:);
n_prev = 2;
if lvl>1
    if lvl==1.5 % almost quantum
        for k=1:n_a
            for l=1:n_b
                T{end+1,1} = k;
                T{end,2} = l;
            end
        end
    elseif lvl>1.5
        for N=2:lvl
            T_prev = T(n_prev:end,:); % previous hierachy
            n_prev = size(T,1)+1;
            for k=1:size(T_prev,1)
                for l=1:size(T0,1)
                    t_new = CombT(T_prev(k,:),T0(l,:),n_xy,n_ab);
                    if isempty(t_new{1})==0 && isempty(t_new{2})==0
                        if DoesExist(T,t_new)<=0 % if it doesn't exist yet
                            T = [T;t_new];
                        end
                    end
                end
            end
        end
        % now remove the same moments taking into account conjugation
        list = 1:size(T,1);
        for k=1:size(T,1)
            list = list(list~=k);
            if DoesExist(T(list,:),T(k,:))==0 % if unique
                list = sort([k,list]);
            end
        end
        T = T(list,:);
        % and let's sort by the size of terms
        for k=1:size(T,1)
            n_sort(k) = numel(T{k,1})+numel(T{k,2});
        end
        [~,list] = sort(n_sort);
        T = T(list,:);
    end                      
end
% now let's constract a matrix
M = zeros(size(T,1));
V_o{1,1} = 0;
V_o{1,2} = 0;
V_o = [V_o;T0]; % observed
for k=1:n_a
    for l=1:n_b
        V_o{end+1,1} = k;
        V_o{end,2} = l;
    end
end
V_u = []; % unobserved
for k=1:size(T,1) % let k count columns
    for l=1:k
        t_kl = CombT(T(k,:), FlipT(T(l,:)),n_xy,n_ab);
        if isempty(t_kl{1}) || isempty(t_kl{2})
            M(l,k) = 0;
        else
            if numel(t_kl{1})+numel(t_kl{2})==2 % observed
                M(l,k) = DoesExist(V_o,t_kl);
            else
                f = DoesExist(V_u,t_kl);
                if f~=0
                    M(l,k) = 1i*f;
                else
                    M(l,k) = 1i*(size(V_u,1)+1);
                    V_u = [V_u; t_kl];
                end
            end      
        end
    end
end
M = M+M'-diag(diag(M'));
% need to fix that diagonal unobserved elements don't have imaginary part
dIm = diag(imag(M));
dIm = dIm(dIm~=0);
for k=1:numel(dIm)
    M(abs(imag(M))==dIm(k))=1i*dIm(k);
end
end


function t = FlipT(T)
t{1} = fliplr(T{1});
t{2} = fliplr(T{2});
end
    
function T = CombT(T1,T2,n_xy,n_ab)
% Combine two moments
% e.g. {[0], [2 1]} and {[2], [2]} gives {[2], [2 1 2]},
% but {[1 2], [2 3]} and {[2 3], [3 4 2]} gives {[1 2 3], [2 3 4 2]}
T{1,1} = 0;
T{1,2} = 0;
for k=1:2
    if nnz(T1{k})==0
        T{k} = T2{k};
    elseif nnz(T2{k})==0
        T{k} = T1{k};
    else
        if isOrt(T1{k}(end),T2{k}(1),n_xy(k),n_ab(k))
            T{k} = [];
        else
            if T1{k}(end)==T2{k}(1)
                T{k} = [T1{k},T2{k}(2:end)];
            else
                T{k} = [T1{k},T2{k}];
            end
        end
    end  
end
end

function e = isOrt(p1,p2,n_s,n_o)
% checks is projectors p1 and p2 are orthogonal
o = i2s([n_s,n_o-1],[p1;p2]);
e = sign(p1*p2)*isequal(o(1,1),o(2,1))*(1-isequal(o(1,2),o(2,2)));
end

function f = DoesExist(T,t)
% finds t in T
f = 0;
for k=1:size(T,1)
    if isequal(T{k,1},t{1}) && isequal(T{k,2},t{2})
        f = k;
    elseif isequal(fliplr(T{k,1}),t{1}) && isequal(fliplr(T{k,2}),t{2})
        f = -k;
    end
end
end