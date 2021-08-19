function [pressao_telhado_x, pressao_telhado_y] = pressao_telhado(pressao, delta)
contador = 0;                                           %contador para a vizinhanca do noh
pressao_telhado = zeros(3/delta+1, 6/delta+1);          %matriz de zeros com tamanho igual ao redor do telhado
pressao_telhado_x = pressao_telhado;                    %matriz para componente horizontal da pressao
pressao_telhado_y = pressao_telhado;                    %matriz para componente vertical da pressao

for j = (15/delta+1):(21/delta+1)                       %loop que percorre psi(i,j) apenas na parte ao redor do telhado
    for i = (3/delta+1):(6/delta+1)
        
        if pressao(i, j) == 0 && (pressao(i+1, j+1)~=0 || pressao(i+1, j-1)~= 0)  %se o ponto for imediatamente acima do telhado pois o ponto lateral ou abaixo esta dentro do telhado                                    
            for j_local = (j-1):(j+1)                   %cria variaveis que checam a vizinhanca do no atual
                for i_local = i:(i+1)
                    
                    if pressao(i_local, j_local) ~= 0   %se a pressao no noh vizinho for diferente de 0
                        pressao_telhado(i-3/delta, j-15/delta) = pressao(i_local, j_local) + pressao_telhado(i-3/delta, j-15/delta); %soma as pressoes entre nos vizinhos
                        contador = contador + 1;        %atualiza contador (para realizar a media)
                        
                        if j == 18/delta + 1            %calcula angulo para dividir a pressao em componentes x e y
                            theta = -pi/2;              
                        else    
                            theta = atan((i -(6/delta+1))/(j - (18/delta+1)))-pi/2;
                        end
                    end
                end
            end
            pressao_telhado(i-3/delta, j-15/delta) = pressao_telhado(i-3/delta, j-15/delta)/contador;       %calcula media entre somatoria das pressoes da vizinhanca e quantidade de nos da vizinhanca nao nulos
            pressao_telhado_x(i-3/delta, j-15/delta) = pressao_telhado(i-3/delta, j-15/delta)*cos(theta);   %calcula componente horizontal da pressao
            pressao_telhado_y(i-3/delta, j-15/delta) = pressao_telhado(i-3/delta, j-15/delta)*sin(theta);   %calcula componente vertical da pressao
            contador = 0;
        end
    end
end
pressao_telhado_x(3/delta, 3/delta+1) = 0;          %deixa o noh imediamente abaixo do topo do telhado com pressoes nulas
pressao_telhado_y(3/delta, 3/delta+1) = 0;          %para nao prejudicar a plotagem da imagem
end



