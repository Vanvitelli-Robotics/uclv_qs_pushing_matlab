
stl  = stlread("../../cad_models/cad_santal_centered1.stl");
P = stl.Points;
z_limit = 0.1;
Pxy = P(abs(P(:,3))<z_limit,:);
Pxy = Pxy(:,1:2);

[~, ind] = min(Pxy(:,1));
Pxy_tmp = Pxy(ind,:);
Pxy(ind,:) = Inf;

Pxy_sorted = zeros(size(Pxy));
Pxy_sorted(1,:) = Pxy_tmp;

for i = 2:length(Pxy)
    [~, ind_tmp] = min(vecnorm((Pxy-Pxy_tmp)'));
    Pxy_tmp = Pxy(ind_tmp,:);
    Pxy_sorted(i,:) = Pxy_tmp;
    Pxy(ind_tmp,:) = Inf;
end

