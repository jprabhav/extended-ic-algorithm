function [Mp,Mr] = SplitMatrix(M)
ro = unique(real(M)); % observed
ro = ro(ro~=0);
ru_r = unique(imag(M(imag(M)>0))); % unobserved real
%ru_i = unique(imag(M(imag(M)<0))); % unobserved imaginary
%ru_i = flipud(ru_i);
Mp = zeros(size(M,1),size(M,2),numel(ro));
Mr = zeros(size(M,1),size(M,2),numel(ru_r));
%Mi = zeros(size(M,1),size(M,2),numel(ru_i));
for k=1:numel(ro)
    Mk = zeros(size(M));
    Mk(real(M)==ro(k)) = 1;
    Mp(:,:,k) = Mk;
end
for k=1:numel(ru_r)
    Mk = zeros(size(M));
    Mk(imag(M)==ru_r(k)) = 1;
    Mr(:,:,k) = sign(Mk+Mk');
end
%for k=1:numel(ru_i)
%    Mk = zeros(size(M));
%    Mk(imag(M)==ru_i(k)) = -1i;
%    Mi(:,:,k) = sign(Mk+Mk');
%end
end