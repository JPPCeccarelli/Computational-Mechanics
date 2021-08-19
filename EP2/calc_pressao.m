function [pressao] = calc_pressao(u, v, delta)
Ni = 24/delta + 1;         %constantes 
Nj = 36/delta +1;
gamma = 1.4;
ro = 1.25;
pressao = zeros(Ni, Nj);   %matriz de pressao em cada noh

    for i = 1:Ni
        for j = 1:Nj
            
            %calcula pressao no dominio em relacao ao telhado, usando a formula do enunciado
            pressao(i, j) = (ro/gamma)*(gamma - 1)*(-(u(i,j)^2 + v(i,j)^2))/2; 
            
        end
    end
end

