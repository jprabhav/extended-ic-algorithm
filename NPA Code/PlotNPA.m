function [alpha_ml,beta_ml,alpha_aq,beta_aq,alpha_qm,beta_qm] = PlotNPA(d,nn)
% COMPAREMLANDIC compares correlated IC vs ML and AQ
% This code needs YALMIP and MOSEK (or other solver, but that needs to be changed in PlotValarea)
% There is no display of progress!
%
% nn -- number of points
n_xy = [d,2];

%q = 1;
%p0 = q*PRbox(d,0)+(1-q)*(ones(1, prod(n_xy)*d^2) / d^2);
p0 = PRbox(d,0);

% local boxes (I guess you want to search over possible configurations here)

p1 = 1/3*(LocalBoxd(n_xy,d,[2,0,0,0])+LocalBoxd(n_xy,d,[2,1,0,2])+LocalBoxd(n_xy,d,[2,2,0,1]));
%p1 =LocalBoxd(n_xy,d,[2,0,0,0]);

p2 = ones(1, prod(n_xy)*d^2) / d^2;
%p2 = 1/3*(LocalBoxd(n_xy,d,[2,0,0,0])+LocalBoxd(n_xy,d,[2,1,0,2])+LocalBoxd(n_xy,d,[2,2,0,1]));

[alpha_ml,beta_ml] = PlotValarea(d,nn,1,p0,p1,p2);
[alpha_aq,beta_aq] = PlotValarea(d,nn,1.5,p0,p1,p2);
[alpha_qm,beta_qm] = PlotValarea(d,nn,2,p0,p1,p2);

end


%{
function p = LocalBoxd(n_xy,d,l)
p = zeros(1,prod(n_xy)*d^2);

for x=0:n_xy(1)-1
    for y=0:n_xy(2)-1
        for a=0:d-1
            for b=0:d-1
                if a==mod(x*l(1)+l(2),d) && b==mod(y*l(3)+l(4),d)
                    for j=0:d-1
                        p(s2i([d d 2 d],[x mod(a+j,d) y mod(b-j,d)]+1)) = 1/d;
                    end
                end
            end
        end
    end
end
end
%}

function p = LocalBoxd(n_xy,d,l)
p = zeros(1,prod(n_xy)*d^2);
for x=0:n_xy(1)-1
    for y=0:n_xy(2)-1
        for a=0:d-1
            for b=0:d-1
                % This condition identifies the unique deterministic outcome (a,b)
                if a==mod(x*l(1)+l(2),d) && b==mod(y*l(3)+l(4),d)
                    % We set the probability to 1 instead of summing over j
                    p(s2i([d d 2 d],[x a y b]+1)) = 1;
                end
            end
        end
    end
end
end

