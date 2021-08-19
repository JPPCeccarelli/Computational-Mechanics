function [forca] = calc_forca(delta, pressaoty)

%area = pi*R*L (m2)
area = pi*3*60;         

%forca inicialmente eh nula para ir somando componentes de pressao
forca = 0;                                 

for i = 1:(3/delta+1)
    for j = 1:(6/delta+1)
        
        %vai somando cada valor de cada ponto com pressao no telhado
        forca = pressaoty(i,j) + forca;  
        
    end
end

%multiplica soma das pressoes por area total
forca = forca*area;                         