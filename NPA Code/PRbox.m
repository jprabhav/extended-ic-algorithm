function p = PRbox(d, s)
% PRBOX generates the Popescu-Rohrlich box for dimension d
%   d   -- local dimension (number of outputs)
%   s   -- shift parameter (0 = standard PR box)
%
% Scenario: Alice has d inputs (x=0..d-1), Bob has 2 inputs (y=0,1)
% Distribution: P(a,b|x,y) = 1/d if a+b = x*y+s (mod d), else 0
% Indexing convention: s2i([d d 2 d], [x+1 a+1 y+1 b+1])

n_xy = [d, 2];
p = zeros(1, prod(n_xy) * d^2);  % length = 2*d^3

for x = 0:d-1
    for y = 0:1
        for a = 0:d-1
            b = mod(x*y + s - a, d);
            p(s2i([d d 2 d], [x+1, a+1, y+1, b+1])) = 1/d;
        end
    end
end
end