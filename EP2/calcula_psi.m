function [psi] = calcula_psi(delta, lambda, converg)
    Ni = 24/delta + 1;              %Divisao da malha em (36/delta+1)x(24/delta+1) pontos
    Nj = 36/delta + 1;              %obs: para gerar graficos mais facilmente invertemos os eixos em relacao as equacoes anteriormente vistas no equacionamento analitico
    psi = zeros(Ni, Nj);            %funcao de corrente psi(i, j) inicialmente nula
    psi_velho = psi;                %psi_velho serve para armazenar dado anterior, para aplicar sobrerrelaxacao
    eps = ones(Ni, Nj);             %epsilon (taxa de erro), inicialmente com 1 e vai diminuindo ateh 0.01
    eps_maximo = 1;                 %valor maximo de epsilon, para quebrar o loop quando for menor que 0.01
    
    while eps_maximo >= converg     %loop roda ate eps_maximo, em todos os pontos da malha, for menor que 0.01
    for i = 1:Ni
        for j = 1:Nj                   
            
           %chão, identificado por todos os pontos com i == 1 e j qualquer
           if i == 1
               eps(i, j) = 0;       %nao faz nada nesse ponto com psi pois o vento nao flui no chao
                                    %mas mesmo assim atualiza o epsilon
                                    %para 0
           end
        
            %ponto na extrema esquerda
           if (j == 1)&&(i == Ni)
                psi_velho(i, j) = psi(i, j);                               %aplica armazena valor anterior
                psi(i, j) = (100*delta/3.6 + psi(i-1, j) + psi(i, j+1))/2; %aplica equacoes deduzidas para o ponto na extrema esquerda
                psi(i, j) = lambda*psi(i, j) + (1-lambda)*psi_velho(i, j); %aplica sobrerrelaxacao
                eps(i, j) = abs((psi(i, j) - psi_velho(i, j))/psi(i, j));  %atualiza o erro neste noh
           end
        
          %a partir deste ponto nao iremos comentar as equacoes analogas a de cima 
          %(borda esquerda, borda direita, topo, ponto extremo direito), pois teriam comentarios praticamente identicos
          
            %borda esquerda
           if(j == 1)&&(i~=1)&&(i~=Ni)
                psi_velho(i, j) = psi(i, j);
                psi(i, j) = (psi(i-1, j)+psi(i+1,j)+2*psi(i,j+1))/4;
                psi(i, j) = lambda*psi(i, j) + (1-lambda)*psi_velho(i, j);
                eps(i, j) = abs((psi(i, j) - psi_velho(i, j))/psi(i, j));
            end
        
             %topo
            if(i == Ni)&&(j ~= 1)&&(j ~= Nj)
                psi_velho(i, j) = psi(i, j);
                psi(i, j) = (psi(i, j-1)+psi(i, j+1)+2*psi(i-1, j)+200*delta/3.6)/4;
                psi(i, j) = lambda*psi(i, j) + (1-lambda)*psi_velho(i, j);
                eps(i, j) = abs((psi(i, j) - psi_velho(i, j))/psi(i, j));   
            end
            
        %ponto extremo a direita
           if(i == Ni)&&(j==Nj)
               psi_velho(i, j) = psi(i, j);
               psi(i, j) = (100*delta/3.6 + psi(i-1, j) + psi(i, j-1))/2;
               psi(i, j) = lambda*psi(i, j) + (1-lambda)*psi_velho(i, j);
               eps(i, j) = abs((psi(i, j) - psi_velho(i, j))/psi(i, j));      
           end
           
        %lado direito
           if(j == Nj)&&(i ~= 1)&&(i~=Ni)
               psi_velho(i, j) = psi(i, j);
               psi(i, j) = ( psi(i+1, j) + psi(i-1, j) + 2*psi(i, j-1))/4;
               psi(i, j) = lambda*psi(i, j) + (1-lambda)*psi_velho(i, j);
               eps(i, j) = abs((psi(i, j) - psi_velho(i, j))/psi(i, j));     
           end
            
           %hangar e resto da malha
            if(j > 1 && j < Nj && i > 1 && i < Ni)          %condicao que verificar se os pontos estao no centro da malha
                %meio da malha
                if(j >= 15/delta+1 && j <= 21/delta+1)      %condicao que verifica se j esta na parte do meio do hangar
                    %parte quadrada do hangar
                    if(i > 1 && i < 3/delta+1)              %parte quadrada do hangar em que nao ha fluxo de vento 
                                                            %e tambem nao ha condicao de borda irregular
                        eps(i, j) = 0;                      %apenas atualiza erro para 0
                        
                    %telhado redondo, pontos imediatamente acima dele
                    elseif(i <= 6/delta+1 && i >= 3/delta+1)
                        %lado direito imediamente superior ao telhado, ptos
                        %ao lado e embaixo tocando
                        
                        if(((i-3/delta-1)^2+(j-2-18/delta)^2 <= (3/delta)^2)&&((i-3/delta-2)^2+(j-1-18/delta)^2 <= (3/delta)^2)&&((i-3/delta-1)^2+(j-1-18/delta)^2 > (3/delta)^2))
                            %condicao que verifica se os pontos (i-1, j) e
                            %(i, j-1) estao dentro do telhado
                                                                                %como visto no equacionamento analitico,
                            theta_y = acos(abs((j - (18/delta +1))*delta/3));   %calcula distancia entre o ponto imediatamente fora do telhado e
                            dif_y = i-(3/delta+1) - 3*sin(theta_y)/delta;       %o telhado, para aplicar cond de contorno irregular em y
                            theta_x = asin((i-(3/delta+1))*delta/3);            %calcula distancia entre o ponto imediatamente fora do telhado e
                            dif_x = abs(j-(18/delta+1)) - 3*cos(theta_x)/delta; %o telhado, para aplicar cond de contorno irregular em y
                            a = (dif_y*delta^2 + dif_y^2*delta);                %a e b sao simplificacoes para diminuir o tamanho das equacoes
                            b = (dif_x*delta^2 + dif_x^2*delta);
                            psi_velho(i, j) = psi(i, j);                        %armazena valor de psi(i, j)
                            psi(i, j) = (b*dif_y*psi(i+1,j)+a*dif_x*psi(i,j+1))/(b*(dif_y+delta)+a*(dif_x+delta));  %aplica equacao deduzida
                            psi(i, j) = lambda*psi(i, j) + (1-lambda)*psi_velho(i, j);  %sobrerrelaxacao
                            eps(i, j) = abs((psi(i, j) - psi_velho(i, j))/psi(i, j));   %verificacao do erro
                            
                        %lateral encosta mas embaixo nao, lado direito
                        elseif(((i-3/delta-1)^2+(j-18/delta-2)^2 <= (3/delta)^2)&&((i-3/delta-1)^2+(j-1-18/delta)^2 > (3/delta)^2)&&((i-3/delta-2)^2+(j-1-18/delta)^2 > (3/delta)^2))
                            %condicao que verifica se apenas o ponto
                            %(i, j-1) esta dentro do telhado
                            theta = asin((i-(3/delta+1))*delta/3);
                            dif = abs(j-(18/delta+1)) - 3*cos(theta)/delta;
                            b = (dif^2+delta*dif);
                            psi_velho(i, j) = psi(i, j);
                            psi(i, j) = (b*(psi(i+1,j)+psi(i-1,j))+2*delta*dif*psi(i,j+1))/(2*((dif+delta)^2));
                            psi(i, j) = lambda*psi(i, j) + (1-lambda)*psi_velho(i, j); 
                            eps(i, j) = abs((psi(i, j) - psi_velho(i, j))/psi(i, j));
                            
                        %lateral nao encosta mas embaixo sim, lado direito
                        elseif(((i-3/delta-1)^2+(j-18/delta-2)^2 > (3/delta)^2)&&((i-3/delta-1)^2+(j-1-18/delta)^2 > (3/delta)^2)&&((i-3/delta-2)^2+(j-1-18/delta)^2 <= (3/delta)^2))
                            %condicao que verifica se apenas o ponto
                            %(i-1, j) esta dentro do telhado
                            theta = acos(abs((j - (18/delta +1))*delta/3));
                            dif = i-(3/delta+1) - 3*sin(theta)/delta;
                            b = (dif^2+delta*dif);
                            psi_velho(i, j) = psi(i, j);
                            psi(i, j) = (b*(psi(i,j+1)+psi(i,j-1))+2*delta*dif*psi(i+1,j))/(2*((dif+delta)^2));
                            psi(i, j) = lambda*psi(i, j) + (1-lambda)*psi_velho(i, j); 
                            eps(i, j) = abs((psi(i, j) - psi_velho(i, j))/psi(i, j));
                            
                        %lado esquerdo, tambem imediatamento superior, ptos
                        %ao lado e embaixo tocando
                        elseif(((i-3/delta-1)^2+(j-18/delta)^2 <= (3/delta)^2)&&((i-3/delta-2)^2+(j-1-18/delta)^2 <= (3/delta)^2)&&((i-3/delta-1)^2+(j-1-18/delta)^2 > (3/delta)^2))
                            %condicao que verifica se os pontos (i-1, j) e
                            %(i, j+1) estao dentro do telhado
                            theta_y = acos(abs((j - (18/delta +1))*delta/3));
                            dif_y = i-(3/delta+1) - 3*sin(theta_y)/delta;
                            theta_x = asin((i-(3/delta+1))*delta/3);
                            dif_x = abs(j-(18/delta+1)) - 3*cos(theta_x)/delta;
                            a = (dif_y*delta^2 + dif_y^2*delta);
                            b = (dif_x*delta^2 + dif_x^2*delta);
                            psi_velho(i, j) = psi(i, j);
                            psi(i, j) = (b*dif_y*psi(i+1,j)+a*dif_x*psi(i,j-1))/(b*(dif_y+delta)+a*(dif_x+delta));
                            psi(i, j) = lambda*psi(i, j) + (1-lambda)*psi_velho(i, j); 
                            eps(i, j) = abs((psi(i, j) - psi_velho(i, j))/psi(i, j));
                            
                        %lateral encosta mas embaixo nao, lado esquerdo
                        elseif(((i-3/delta-1)^2+(j-18/delta)^2 <= (3/delta)^2)&&((i-3/delta-1)^2+(j-1-18/delta)^2 > (3/delta)^2)&&((i-3/delta-2)^2+(j-1-18/delta)^2 > (3/delta)^2))
                            %condicao que verifica se apenas o ponto
                            %(i, j+1) esta dentro do telhado
                            theta = asin((i-(3/delta+1))*delta/3);
                            dif = abs(j-(18/delta+1)) - 3*cos(theta)/delta;
                            b = (dif^2+delta*dif);
                            psi_velho(i, j) = psi(i, j);
                            psi(i, j) = (b*(psi(i+1,j)+psi(i-1,j))+2*delta*dif*psi(i,j-1))/(2*((dif+delta)^2));
                            psi(i, j) = lambda*psi(i, j) + (1-lambda)*psi_velho(i, j); 
                            eps(i, j) = abs((psi(i, j) - psi_velho(i, j))/psi(i, j));
                            
                        %lateral nao encosta mas embaixo sim, lado esquerdo
                        elseif(((i-3/delta-1)^2+(j-18/delta)^2 > (3/delta)^2)&&((i-3/delta-1)^2+(j-1-18/delta)^2 > (3/delta)^2)&&((i-3/delta-2)^2+(j-1-18/delta)^2 <= (3/delta)^2))
                            %condicao que verifica se apenas o ponto
                            %(i-1, j) esta dentro do telhado
                            theta = acos(abs((j - (18/delta +1))*delta/3));
                            dif = i-(3/delta+1) - 3*sin(theta)/delta;
                            b = (dif^2+delta*dif);
                            psi_velho(i, j) = psi(i, j);
                            psi(i, j) = (b*(psi(i,j+1)+psi(i,j-1))+2*delta*dif*psi(i+1,j))/(2*((dif+delta)^2));
                            psi(i, j) = lambda*psi(i, j) + (1-lambda)*psi_velho(i, j); 
                            eps(i, j) = abs((psi(i, j) - psi_velho(i, j))/psi(i, j));
                            
                        %interior do telhado, pontos abaixo da
                        %superficie
                        elseif((i-3/delta-1)^2+(j-1-18/delta)^2 <= (3/delta)^2)
                            %atualiza erro como nulo pois nao tem fluxo de ar
                            eps(i, j) = 0;              
                        
                        %os pontos abaixo sao todos para o restante da malha em que nao ha condicoes de contorno
                        %portanto, pode-se aplicar as condicoes de
                        %diferencas finitas centrais
                            
                        %restante acima do telhado, entre 3/8 e 5/8 da
                        %malha
                        else
                            psi_velho(i, j) = psi(i, j);
                            psi(i, j) = ( psi(i+1, j) + psi(i-1, j) + psi(i, j-1)+psi(i, j+1))/4;   %equacao de dif finitas sem condicao de contorno
                            psi(i, j) = lambda*psi(i, j) + (1-lambda)*psi_velho(i, j); 
                            eps(i, j) = abs((psi(i, j) - psi_velho(i, j))/psi(i, j));
                        end
                        
                    %Restante acima do telhado e acima de 5/8 da malha    
                    else
                        psi_velho(i, j) = psi(i, j);
                        psi(i, j) = ( psi(i+1, j) + psi(i-1, j) + psi(i, j-1)+psi(i, j+1))/4;   
                        psi(i, j) = lambda*psi(i, j) + (1-lambda)*psi_velho(i, j); 
                        eps(i, j) = abs((psi(i, j) - psi_velho(i, j))/psi(i, j));
                    end
                    
                %Restante da malha, laterais ao hangar
                else
                    psi_velho(i, j) = psi(i, j);
                    psi(i, j) = ( psi(i+1, j) + psi(i-1, j) + psi(i, j-1)+psi(i, j+1))/4;   
                    psi(i, j) = lambda*psi(i, j) + (1-lambda)*psi_velho(i, j); 
                    eps(i, j) = abs((psi(i, j) - psi_velho(i, j))/psi(i, j));
                end
            end
        end
    end
    
    %atualiza o valor maximo do erro epsilon para a matriz de psi
    eps_maximo = max(eps, [], 'all');          
    end 
end

