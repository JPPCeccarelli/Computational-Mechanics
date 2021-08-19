function [taxa_calor] = calc_taxa_calor(delta,fluxo)
area = pi*3*60; %calcula area do telhado
taxa_calor = 0;
Ni = 24/delta + 1;
Nj = 36/delta + 1;

for i = 1:Ni
    for j = 1:Nj
            taxa_calor = fluxo(i,j)*area + taxa_calor;
        end
    end
end