function [u, v] = velocidades(psi, delta)
Ni = 24/delta + 1;
Nj = 36/delta + 1;
u = zeros(Ni, Nj);  %matriz de velocidades u em cada ponto
v = u;              %matriz de velocidades v em cada ponto

    for i = 1:Ni
        for j = 1:Nj                   

            %chão: nao faz nada pois as velocidades sao 0
            if i == 1

            end

            %pontos extremos
            if ((j == 1)&&(i == Ni))||(j == Nj && i == Ni)  
                %aplica condicao de contorno de borda (v = 0) e topo
                %(u=100/3.6)
                v(i, j) = 0;           
                u(i, j) = 100/3.6;
            end

            %bordas
            if((j == 1)&&(i~=1)&&(i~=Ni))||((j == Nj)&&(i~=1)&&(i~=Ni))
                %aplica condicao de contorno de borda v = 0
                v(i, j) = 0;        
                %aplica eq de primeira diferenca central
                u(i, j) = (psi(i+1, j) - psi(i-1, j))/(2*delta);    
            end

             %topo
            if(i == Ni)&&(j ~= 1)&&(j ~= Nj)
                 %condicao de contorno u = V
                u(i, j) = 100/3.6;       
                %aplica eq de primeira diferenca central
                v(i, j) = (psi(i, j - 1) - psi(i, j + 1))/(2*delta); 
            end
            
            %abaixo temos a mesma ideia: se um dos pontos (lateral ou
            %inferior) encostarem entao aplica-se a serie de taylor 
            %correspondente de ordem 2 tal como explicitada no 
            %equacionamento analitico. Caso nao encoste entao eh aplicada
            %a primeira diferenca central 
            
           %hangar e resto da malha
            if(j > 1 && j < Nj && i > 1 && i < Ni)
                %meio da malha
                if(j >= 15/delta+1 && j <= 21/delta+1)
                    %parte quadrada do hangar
                    if(i > 1 && i < 3/delta+1)
                       
                    %telhado redondo
                    elseif(i <= 6/delta+1 && i >= 3/delta+1)
                        
                        %lado direito imediamente superior ao telhado,
                        %ambos encostam
                        if(((i-3/delta-1)^2+(j-2-18/delta)^2 <= (3/delta)^2)&&((i-3/delta-2)^2+(j-1-18/delta)^2 <= (3/delta)^2)&&((i-3/delta-1)^2+(j-1-18/delta)^2 > (3/delta)^2))
                            theta_y = acos(abs((j - (18/delta +1))*delta/3));
                            dif_y = i-(3/delta+1) - 3*sin(theta_y)/delta;
                            theta_x = asin((i-(3/delta+1))*delta/3);
                            dif_x = abs(j-(18/delta+1)) - 3*cos(theta_x)/delta;
                            u(i,j) = (dif_y^2*psi(i+1,j)+(delta^2-dif_y^2)*psi(i,j))/(dif_y*(delta^2+dif_y*delta));
                            v(i,j) = -(dif_x^2*psi(i,j+1)+(delta^2-dif_x^2)*psi(i,j))/(dif_x*(delta^2+dif_x*delta));
                            
                        %lateral encosta mas embaixo nao, lado direito
                        elseif(((i-3/delta-1)^2+(j-18/delta-2)^2 <= (3/delta)^2)&&((i-3/delta-1)^2+(j-1-18/delta)^2 > (3/delta)^2)&&((i-3/delta-2)^2+(j-1-18/delta)^2 > (3/delta)^2))
                            theta_x = asin((i-(3/delta+1))*delta/3);
                            dif_x = abs(j-(18/delta+1)) - 3*cos(theta_x)/delta;
                            v(i,j) = -(dif_x^2*psi(i,j+1)+(delta^2-dif_x^2)*psi(i,j))/(dif_x*(delta^2+dif_x*delta));
                            u(i,j) = (psi(i+1,j)-psi(i-1,j))/(2*delta);
                            
                        %lateral nao encosta mas embaixo sim, lado direito
                        elseif(((i-3/delta-1)^2+(j-18/delta-2)^2 > (3/delta)^2)&&((i-3/delta-1)^2+(j-1-18/delta)^2 > (3/delta)^2)&&((i-3/delta-2)^2+(j-1-18/delta)^2 <= (3/delta)^2))
                            theta_y = acos(abs((j - (18/delta +1))*delta/3));
                            dif_y = i-(3/delta+1) - 3*sin(theta_y)/delta;
                            v(i,j) = (psi(i,j-1)-psi(i,j+1))/(2*delta);
                            u(i,j) = (dif_y^2*psi(i+1,j)+(delta^2-dif_y^2)*psi(i,j))/(dif_y*(delta^2+dif_y*delta));
                            
                        %lado esquerdo, tambem imediatamento superior,
                        %ambos encostam
                        elseif(((i-3/delta-1)^2+(j-18/delta)^2 <= (3/delta)^2)&&((i-3/delta-2)^2+(j-1-18/delta)^2 <= (3/delta)^2)&&((i-3/delta-1)^2+(j-1-18/delta)^2 > (3/delta)^2))
                            theta_y = acos(abs((j - (18/delta +1))*delta/3));
                            dif_y = i-(3/delta+1) - 3*sin(theta_y)/delta;
                            theta_x = asin((i-(3/delta+1))*delta/3);
                            dif_x = abs(j-(18/delta+1)) - 3*cos(theta_x)/delta;
                            u(i,j) = (dif_y^2*psi(i+1,j)+(delta^2-dif_y^2)*psi(i,j))/(dif_y*(delta^2+dif_y*delta));
                            v(i,j) = (dif_x^2*psi(i,j-1)+(delta^2-dif_x^2)*psi(i,j))/(dif_x*(delta^2+dif_x*delta));
                            
                        %lateral encosta mas embaixo nao, lado esquerdo
                        elseif(((i-3/delta-1)^2+(j-18/delta)^2 <= (3/delta)^2)&&((i-3/delta-1)^2+(j-1-18/delta)^2 > (3/delta)^2)&&((i-3/delta-2)^2+(j-1-18/delta)^2 > (3/delta)^2))
                            theta_x = asin((i-(3/delta+1))*delta/3);
                            dif_x = abs(j-(18/delta+1)) - 3*cos(theta_x)/delta;
                            v(i,j) = (dif_x^2*psi(i,j-1)+(delta^2-dif_x^2)*psi(i,j))/(dif_x*(delta^2+dif_x*delta));
                            u(i,j) = (psi(i+1,j)-psi(i-1,j))/(2*delta);
                            
                        %lateral nao encosta mas embaixo sim, lado esquerdo
                        elseif(((i-3/delta-1)^2+(j-18/delta)^2 > (3/delta)^2)&&((i-3/delta-1)^2+(j-1-18/delta)^2 > (3/delta)^2)&&((i-3/delta-2)^2+(j-1-18/delta)^2 <= (3/delta)^2))
                            theta_y = acos(abs((j - (18/delta +1))*delta/3));
                            dif_y = i-(3/delta+1) - 3*sin(theta_y)/delta;
                            v(i,j) = (psi(i,j-1)-psi(i,j+1))/(2*delta);
                            u(i,j) = (dif_y^2*psi(i+1,j)+(delta^2-dif_y^2)*psi(i,j))/(dif_y*(delta^2+dif_y*delta));
                            
                        %interior do telhado
                        elseif((i-3/delta-1)^2+(j-1-18/delta)^2 <= (3/delta)^2)

                        %restante acima do telhado, entre 3/8 e 5/8 da
                        %malha
                        else
                    		u(i, j) = (psi(i + 1, j) - psi(i -1, j))/(2*delta);
                    		v(i, j) = (psi(i, j-1) - psi(i, j+1))/(2*delta);
                        end
                        
                    %Restante acima do telhado e acima de 5/8 da malha    
                    else
                        u(i, j) = (psi(i + 1, j) - psi(i -1, j))/(2*delta);
                        v(i, j) = (psi(i, j-1) - psi(i, j+1))/(2*delta);
                    end
                    
                %Restante da malha, laterais ao hangar
                else
                    u(i, j) = (psi(i +1, j) - psi(i -1, j))/(2*delta);
                    v(i, j) = (psi(i, j-1) - psi(i, j+1))/(2*delta);
                end
            end
        end
    end
end        
        
