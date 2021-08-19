function[y_x, y_y, y_m, lambda] = modos_de_vibracao(K, M)
%condições de contorno (nohs de base)
M(:,1) = [];
M(1,:) = [];
M(:,1) = [];
M(1,:) = [];
M(:,3) = [];
M(3,:) = [];

K(:,1) = [];
K(1,:) = [];
K(:,1) = [];
K(1,:) = [];
K(:,3) = [];
K(3,:) = [];

[y, lambda] = eig(K/M);

%frequencias de ressonancia
lambda = sort(diag(sqrt(lambda/(4*pi^2))));

%deslocamentos em cada direção para cada noh
y_x = [zeros(1,39); y(2,:); y(4,:); y(7,:); y(10,:); y(13,:); y(16,:); y(19,:); y(22,:); y(25,:); y(28,:); y(31,:); y(34,:); y(37,:)];
y_y = [zeros(1,39); zeros(1,39); y(5,:); y(8,:); y(11,:); y(14,:); y(17,:); y(20,:); y(23,:); y(26,:); y(29,:); y(32,:); y(35,:); y(38,:)];
y_m = [y(1,:); y(3,:); y(6,:); y(9,:); y(12,:); y(15,:); y(18,:); y(21,:); y(24,:); y(27,:); y(30,:); y(33,:); y(36,:); y(39,:)];
end

