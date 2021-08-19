function [fluxo] = calcula_fluxo(delta, temp, u, v)
Ni = 24/delta + 1;
Nj = 36/delta + 1;
fluxo = zeros(Ni, Nj);
k = 0.026;

for j = 1:Nj
    for i = 1:Ni
        %hangar e resto da malha
        if(j > 1 && j < Nj && i > 1 && i < Ni)
            %meio da malha
            if(j >= 15/delta+1 && j <= 21/delta+1)
                %parte quadrada do hangar
                if(i > 1 && i < 3/delta+1)
                    if(u(i,j)<0 &&v(i,j)<0)
                            fluxo(i,j) = -k*((temp(i,j+1)-temp(i,j))/delta + (temp(i+1,j)-temp(i,j))/delta);
                    end

                    if(u(i,j)>0 &&v(i,j)<0)
                        fluxo(i,j) = -k*((temp(i,j)-temp(i,j-1))/delta + (temp(i+1,j)-temp(i,j))/delta);
                    end

                    if(u(i,j)<0 &&v(i,j)>0)
                        fluxo(i,j) = -k*((temp(i,j+1)-temp(i,j))/delta + (temp(i,j)-temp(i-1,j))/delta);
                    end

                    if(u(i,j)>0 &&v(i,j)>0)
                        fluxo(i,j) = -k*((temp(i,j)-temp(i,j-1))/delta + (temp(i,j)-temp(i-1,j))/delta);
                    end
                %telhado redondo
                elseif(i <= 6/delta+1 && i >= 3/delta+1)
                    %lado direito imediamente superior ao telhado, ptos
                    %ao lado e embaixo tocando 5
                    if(((i-3/delta-1)^2+(j-2-18/delta)^2 <= (3/delta)^2)&&((i-3/delta-2)^2+(j-1-18/delta)^2 <= (3/delta)^2)&&((i-3/delta-1)^2+(j-1-18/delta)^2 > (3/delta)^2))
                        theta_y = acos(abs((j - (18/delta +1))*delta/3));
                        dif_y = i-(3/delta+1) - 3*sin(theta_y)/delta;
                        theta_x = asin((i-(3/delta+1))*delta/3);
                        dif_x = abs(j-(18/delta+1)) - 3*cos(theta_x)/delta;
                        if(u(i,j)<0 &&v(i,j)<0)
                            fluxo(i,j) = -k*((temp(i,j+1)-temp(i,j))*cos(theta_x)/delta + (temp(i+1,j)-temp(i,j))*sin(theta_y)/delta);
                        end

                        if(u(i,j)>0 &&v(i,j)<0)
                            fluxo(i,j) = -k*((temp(i,j)-temp(i,j-1))*cos(theta_x)/dif_x + (temp(i+1,j)-temp(i,j))*sin(theta_y)/delta);
                        end

                        if(u(i,j)<0 &&v(i,j)>0)
                            fluxo(i,j) = -k*((temp(i,j+1)-temp(i,j))*cos(theta_x)/delta + (temp(i,j)-temp(i-1,j))*sin(theta_y)/dif_y);
                        end

                        if(u(i,j)>0 &&v(i,j)>0)
                            fluxo(i,j) = -k*((temp(i,j)-temp(i,j-1))*cos(theta_x)/dif_x + (temp(i,j)-temp(i-1,j))*sin(theta_y)/dif_y);
                        end

                    %lateral encosta mas embaixo nao, lado direito
                    %6
                    elseif(((i-3/delta-1)^2+(j-18/delta-2)^2 <= (3/delta)^2)&&((i-3/delta-1)^2+(j-1-18/delta)^2 > (3/delta)^2)&&((i-3/delta-2)^2+(j-1-18/delta)^2 > (3/delta)^2))
                        theta = asin((i-(3/delta+1))*delta/3);
                        dif_x = abs(j-(18/delta+1)) - 3*cos(theta)/delta;
                        if(u(i,j)<0 &&v(i,j)<0)
                            fluxo(i,j) = -k*((temp(i,j+1)-temp(i,j))*cos(theta)/delta + (temp(i+1,j)-temp(i,j))*sin(theta)/delta);
                        end

                        if(u(i,j)>0 &&v(i,j)<0)
                            fluxo(i,j) = -k*((temp(i,j)-temp(i,j-1))*cos(theta)/dif_x + (temp(i+1,j)-temp(i,j))*sin(theta)/delta);
                        end

                        if(u(i,j)<0 &&v(i,j)>0)
                            fluxo(i,j) = -k*((temp(i,j+1)-temp(i,j))*cos(theta)/delta + (temp(i,j)-temp(i-1,j))*sin(theta)/delta);
                        end

                        if(u(i,j)>0 &&v(i,j)>0)
                            fluxo(i,j) = -k*((temp(i,j)-temp(i,j-1))*cos(theta)/dif_x + (temp(i,j)-temp(i-1,j))*sin(theta)/delta);
                        end
                    %lateral nao encosta mas embaixo sim, lado
                    %direito 4
                    elseif(((i-3/delta-1)^2+(j-18/delta-2)^2 > (3/delta)^2)&&((i-3/delta-1)^2+(j-1-18/delta)^2 > (3/delta)^2)&&((i-3/delta-2)^2+(j-1-18/delta)^2 <= (3/delta)^2))
                        theta = acos(abs((j - (18/delta +1))*delta/3));
                        dif_y = i-(3/delta+1) - 3*sin(theta)/delta;
                        if(u(i,j)<0 &&v(i,j)<0)
                            fluxo(i,j) = -k*((temp(i,j+1)-temp(i,j))*cos(theta)/delta + (temp(i+1,j)-temp(i,j))*sin(theta)/delta);
                        end

                        if(u(i,j)>0 &&v(i,j)<0)
                            fluxo(i,j) = -k*((temp(i,j)-temp(i,j-1))*cos(theta)/delta + (temp(i+1,j)-temp(i,j))*sin(theta)/delta);
                        end

                        if(u(i,j)<0 &&v(i,j)>0)
                            fluxo(i,j) = -k*((temp(i,j+1)-temp(i,j))*cos(theta)/delta + (temp(i,j)-temp(i-1,j))*sin(theta)/dif_y);
                        end

                        if(u(i,j)>0 &&v(i,j)>0)
                            fluxo(i,j) = -k*((temp(i,j)-temp(i,j-1))*cos(theta)/delta + (temp(i,j)-temp(i-1,j))*sin(theta)/dif_y);
                        end

                    %lado esquerdo, tambem imediatamento superior, ptos
                    %ao lado e embaixo tocando 2
                    elseif(((i-3/delta-1)^2+(j-18/delta)^2 <= (3/delta)^2)&&((i-3/delta-2)^2+(j-1-18/delta)^2 <= (3/delta)^2)&&((i-3/delta-1)^2+(j-1-18/delta)^2 > (3/delta)^2))
                        theta_y = acos(abs((j - (18/delta +1))*delta/3));
                        dif_y = i-(3/delta+1) - 3*sin(theta_y)/delta;
                        theta_x = asin((i-(3/delta+1))*delta/3);
                        dif_x = abs(j-(18/delta+1)) - 3*cos(theta_x)/delta;
                        if(u(i,j)<0 &&v(i,j)<0)
                            fluxo(i,j) = -k*((temp(i,j+1)-temp(i,j))*cos(theta_x)/dif_x + (temp(i+1,j)-temp(i,j))*sin(theta_y)/delta);
                        end

                        if(u(i,j)>0 &&v(i,j)<0)
                            fluxo(i,j) = -k*((temp(i,j)-temp(i,j-1))*cos(theta_x)/delta + (temp(i+1,j)-temp(i,j))*sin(theta_y)/delta);
                        end

                        if(u(i,j)<0 &&v(i,j)>0)
                            fluxo(i,j) = -k*((temp(i,j+1)-temp(i,j))*cos(theta_x)/dif_x + (temp(i,j)-temp(i-1,j))*sin(theta_y)/dif_y);
                        end

                        if(u(i,j)>0 &&v(i,j)>0)
                            fluxo(i,j) = -k*((temp(i,j)-temp(i,j-1))*cos(theta_x)/delta + (temp(i,j)-temp(i-1,j))*sin(theta_y)/dif_y);
                        end
                    %lateral encosta mas embaixo nao, lado esquerdo
                    %1
                    elseif(((i-3/delta-1)^2+(j-18/delta)^2 <= (3/delta)^2)&&((i-3/delta-1)^2+(j-1-18/delta)^2 > (3/delta)^2)&&((i-3/delta-2)^2+(j-1-18/delta)^2 > (3/delta)^2))
                        theta = asin((i-(3/delta+1))*delta/3);
                        dif_x = abs(j-(18/delta+1)) - 3*cos(theta)/delta;
                        if(u(i,j)<0 &&v(i,j)<0)
                            fluxo(i,j) = -k*((temp(i,j+1)-temp(i,j))*cos(theta)/delta + (temp(i+1,j)-temp(i,j))*sin(theta)/delta);
                        end

                        if(u(i,j)>0 &&v(i,j)<0)
                            fluxo(i,j) = -k*((temp(i,j)-temp(i,j-1))*cos(theta)/dif_x + (temp(i+1,j)-temp(i,j))*sin(theta)/delta);
                        end

                        if(u(i,j)<0 &&v(i,j)>0)
                            fluxo(i,j) = -k*((temp(i,j+1)-temp(i,j))*cos(theta)/delta + (temp(i,j)-temp(i-1,j))*sin(theta)/delta);
                        end

                        if(u(i,j)>0 &&v(i,j)>0)
                            fluxo(i,j) = -k*((temp(i,j)-temp(i,j-1))*cos(theta)/dif_x + (temp(i,j)-temp(i-1,j))*sin(theta)/delta);
                        end

                    %lateral nao encosta mas embaixo sim, lado
                    %esquerdo 3
                    elseif(((i-3/delta-1)^2+(j-18/delta)^2 > (3/delta)^2)&&((i-3/delta-1)^2+(j-1-18/delta)^2 > (3/delta)^2)&&((i-3/delta-2)^2+(j-1-18/delta)^2 <= (3/delta)^2))
                        theta = acos(abs((j - (18/delta +1))*delta/3));
                        dif_y = i-(3/delta+1) - 3*sin(theta)/delta;
                        if(u(i,j)<0 &&v(i,j)<0)
                            fluxo(i,j) = -k*((temp(i,j+1)-temp(i,j))*cos(theta)/delta + (temp(i+1,j)-temp(i,j))*sin(theta)/delta);
                        end

                        if(u(i,j)>0 &&v(i,j)<0)
                            fluxo(i,j) = -k*((temp(i,j)-temp(i,j-1))*cos(theta)/delta + (temp(i+1,j)-temp(i,j))*sin(theta)/delta);
                        end

                        if(u(i,j)<0 &&v(i,j)>0)
                            fluxo(i,j) = -k*((temp(i,j+1)-temp(i,j))*cos(theta)/delta + (temp(i,j)-temp(i-1,j))*sin(theta)/dif_y);
                        end

                        if(u(i,j)>0 &&v(i,j)>0)
                            fluxo(i,j) = -k*((temp(i,j)-temp(i,j-1))*cos(theta)/delta + (temp(i,j)-temp(i-1,j))*sin(theta)/dif_y);
                        end
                    end
                end
            end
        end
    end
end
end

