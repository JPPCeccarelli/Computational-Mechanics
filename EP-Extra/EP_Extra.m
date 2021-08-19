%Exercicio programa extra
%Autor: Joao P P Ceccarelli / Henrique Yamamoto
%Para alterar a discretizacao basta mudar h na linha 17

clear all;

%dados do enunciado
vmax = 500;
imax = 200;
rc = 0.02;
mu_ar = 1.2566*10^(-6);
mu_solo = 2*1.2567*10^(-6);
sigma_ar = 10^-10;
sigma_solo = 10^-2;

%discretizacao (mudar aqui caso quiser comparar)
h = 1;

%raio e tamanho das matrizes
tam = 52/h+1;
r = 26/h;

%matrizes para guardar nohs dos elementos
%elementos1 = triangulos com maior face virada para cima
%elementos2 = triangulos com maior face virada para baixo
elementos1 = zeros(tam^2);
elementos2 = zeros(tam^2);

%Montagem da matriz de rigidez
K = zeros(tam^2);
Ke = [];
cnt = 1;
for m = 1:tam-1
    for n = 1:tam-1
            %aloca posicoes dos nos
            elementos1(m, n, 1) = n;
            elementos1(m, n, 2) = m;
            elementos1(m, n, 3) = n+h;
            elementos1(m, n, 4) = m;
            elementos1(m, n, 5) = n;
            elementos1(m, n, 6) = m+h;
            
            elementos2(m, n, 1) = n+h;
            elementos2(m, n, 2) = m+h;
            elementos2(m, n, 3) = n;
            elementos2(m, n, 4) = m+h;
            elementos2(m, n, 5) = n+h;
            elementos2(m, n, 6) = m;
            
            %constantes
            b11 = elementos1(m, n, 4) - elementos1(m, n, 6);
            b21 = elementos1(m, n, 6) - elementos1(m, n, 2);
            b31 = elementos1(m, n, 2) - elementos1(m, n, 4); 
            c11 = elementos1(m, n, 5) - elementos1(m, n, 3); 
            c21 = elementos1(m, n, 1) - elementos1(m, n, 5); 
            c31 = elementos1(m, n, 3) - elementos1(m, n, 1); 
            b12 = elementos2(m, n, 4) - elementos2(m, n, 6); 
            b22 = elementos2(m, n, 6) - elementos2(m, n, 2); 
            b32 = elementos2(m, n, 2) - elementos2(m, n, 4); 
            c12 = elementos2(m, n, 5) - elementos2(m, n, 3); 
            c22 = elementos2(m, n, 1) - elementos2(m, n, 5); 
            c32 = elementos2(m, n, 3) - elementos2(m, n, 1);
            
            Ae = (b11*c21 - b21*c11)/2;
        
            Ke = 1/(4*Ae)*[b11*b11+c11*c11 b21*b11+c21*c11 b31*b11+c31*c11;
            b21*b11+c21*c11 b21*b21+c21*c21 b31*b21+c31*c21;
            b31*b11+c31*c11 b21*b31+c21*c31 b31*b31+c31*c31];
        
            K(cnt,cnt) = K(cnt,cnt)+Ke(1,1);
            K(cnt,cnt+1) = K(cnt,cnt+1)+Ke(1,2);
            K(cnt,cnt+tam) = K(cnt,cnt+tam)+Ke(1,3);
            K(cnt+1,cnt) = K(cnt+1,cnt)+Ke(2,1);
            K(cnt+1,cnt+1) = K(cnt+1,cnt+1)+Ke(2,2);
            K(cnt+1,cnt+tam) = K(cnt+1,cnt+tam)+Ke(2,3);
            K(cnt+tam,cnt) = K(cnt+tam,cnt)+Ke(3,1);
            K(cnt+tam,cnt+1) =  K(cnt+tam,cnt+1)+Ke(3,2);
            K(cnt+tam,cnt+tam) = K(cnt+tam,cnt+tam)+Ke(3,3);
            
            Ae = (b12*c22 - b22*c12)/2;
        
            Ke = 1/(4*Ae)*[b12*b12+c12*c12 b22*b12+c22*c12 b32*b12+c32*c12;
            b22*b12+c22*c12 b22*b22+c22*c22 b32*b22+c32*c22;
            b32*b12+c32*c12 b22*b32+c22*c32 b32*b32+c32*c32];
        
            K(cnt+tam+1,cnt+tam+1) = K(cnt+tam+1,cnt+tam+1)+Ke(1,1);
            K(cnt+tam+1,cnt+tam) = K(cnt+tam+1,cnt+tam)+Ke(1,2);
            K(cnt+tam+1,cnt+1) = K(cnt+tam+1,cnt+1)+Ke(1,3);
            K(cnt+tam,cnt+tam+1) = K(cnt+tam,cnt+tam+1)+Ke(2,1);
            K(cnt+tam,cnt+tam) = K(cnt+tam,cnt+tam)+Ke(2,2);
            K(cnt+tam,cnt+1) = K(cnt+tam,cnt+1)+Ke(2,3);
            K(cnt+1,cnt+tam+1) = K(cnt+1,cnt+tam+1)+Ke(3,1);
            K(cnt+1,cnt+tam) = K(cnt+1,cnt+tam)+Ke(3,2);
            K(cnt+1,cnt+1) = K(cnt+1,cnt+1)+Ke(3,3);
            
            cnt = cnt+1;
    end
end

%Condicoes de contorno para potencial eletrico nos pontos conhecidos
%res eh a matriz de carregamentos + CC
%pot_aux armazena a potencia dos nohs num array (tam^2, 1)
%cnt = contador para percorrer tal array
res = zeros(tam^2,1);
pot_aux = ones(tam^2,1);
cnt = 1;
for m = 1:tam
    for n = 1:tam
        if ((m-(tam+1)/2)^2+(n-(tam+1)/2)^2 >= r^2)
            pot_aux(cnt) = 0;
        elseif(m==(ceil(tam/2)+10/h)&&(n==(ceil(tam/2)-6/h)||(n==(ceil(tam/2)+6/h))))||...
                (m==(ceil(tam/2)+14/h)&&(n==(ceil(tam/2)-4/h)||(n==(ceil(tam/2)+4/h))))
            pot_aux(cnt) = vmax;
        end
        cnt = cnt + 1;
    end
end

%Usa algoritmo de eliminacao gaussiana para ir alterando res
%e a matriz de rigidez
for m = 1:tam^2
    for n = 1:tam^2
        %1 eh o valor default da matriz, eh como se nao existisse noh caso
        %o valor seja 1
        if (pot_aux(n) ~= 1)&& m == n
            for k = 1:tam^2
                res(k) = res(k) - K(k,n)*pot_aux(n);
            end
            K(m,:) = 0;
            K(:,n) = 0;
            K(m,n) = 1;
            res(m) = pot_aux(m);
        end
    end
end

%encontra o potencial eletrico (vetor (tam^2, 1))
potencial = linsolve(K,res);

%passa o potencial para o dominio 2D
cnt = 1;
V = zeros(tam);
for m = 1:tam
    for n = 1:tam
        V(m,n) = potencial(cnt);
        cnt = cnt + 1;
    end
end

%condicoes de contorno para o carro
for m = ceil(tam/2):(ceil(tam/2)+2/h)
    for n = (ceil(tam/2)-2/h):(ceil(tam/2)+2/h)
        if(m >= ceil(tam/2) && m <= ceil(tam/2)+1.5/h)&&(n <= ceil(tam/2)+2/h && n>=ceil(tam/2)-2/h)
            V(m,n) = V(ceil(tam/2), ceil(tam/2));
        end
    end
end

cnt = 1;
%matriz de rigidez para variaveis magneticas
%reseta matriz de rigidez antiga 
K = zeros(tam^2);
for m = 1:tam-1
    for n = 1:tam-1            
            b11 = elementos1(m, n, 4) - elementos1(m, n, 6);
            b21 = elementos1(m, n, 6) - elementos1(m, n, 2);
            b31 = elementos1(m, n, 2) - elementos1(m, n, 4); 
            c11 = elementos1(m, n, 5) - elementos1(m, n, 3); 
            c21 = elementos1(m, n, 1) - elementos1(m, n, 5); 
            c31 = elementos1(m, n, 3) - elementos1(m, n, 1); 
            b12 = elementos2(m, n, 4) - elementos2(m, n, 6); 
            b22 = elementos2(m, n, 6) - elementos2(m, n, 2); 
            b32 = elementos2(m, n, 2) - elementos2(m, n, 4); 
            c12 = elementos2(m, n, 5) - elementos2(m, n, 3); 
            c22 = elementos2(m, n, 1) - elementos2(m, n, 5); 
            c32 = elementos2(m, n, 3) - elementos2(m, n, 1); 
            
            Ae = (b11*c21 - b21*c11)/2;
        
            Ke = 1/(4*Ae)*[b11*b11+c11*c11 b21*b11+c21*c11 b31*b11+c31*c11;
            b21*b11+c21*c11 b21*b21+c21*c21 b31*b21+c31*c21;
            b31*b11+c31*c11 b21*b31+c21*c31 b31*b31+c31*c31];
        
            K(cnt,cnt) = K(cnt,cnt)+Ke(1,1);
            K(cnt,cnt+1) = K(cnt,cnt+1)+Ke(1,2);
            K(cnt,cnt+tam) = K(cnt,cnt+tam)+Ke(1,3);
            K(cnt+1,cnt) = K(cnt+1,cnt)+Ke(2,1);
            K(cnt+1,cnt+1) = K(cnt+1,cnt+1)+Ke(2,2);
            K(cnt+1,cnt+tam) = K(cnt+1,cnt+tam)+Ke(2,3);
            K(cnt+tam,cnt) = K(cnt+tam,cnt)+Ke(3,1);
            K(cnt+tam,cnt+1) =  K(cnt+tam,cnt+1)+Ke(3,2);
            K(cnt+tam,cnt+tam) = K(cnt+tam,cnt+tam)+Ke(3,3);
            
            Ae = (b12*c22 - b22*c12)/2;
        
            Ke = 1/(4*Ae)*[b12*b12+c12*c12 b22*b12+c22*c12 b32*b12+c32*c12;
            b22*b12+c22*c12 b22*b22+c22*c22 b32*b22+c32*c22;
            b32*b12+c32*c12 b22*b32+c22*c32 b32*b32+c32*c32];
        
            K(cnt+tam+1,cnt+tam+1) = K(cnt+tam+1,cnt+tam+1)+Ke(1,1);
            K(cnt+tam+1,cnt+tam) = K(cnt+tam+1,cnt+tam)+Ke(1,2);
            K(cnt+tam+1,cnt+1) = K(cnt+tam+1,cnt+1)+Ke(1,3);
            K(cnt+tam,cnt+tam+1) = K(cnt+tam,cnt+tam+1)+Ke(2,1);
            K(cnt+tam,cnt+tam) = K(cnt+tam,cnt+tam)+Ke(2,2);
            K(cnt+tam,cnt+1) = K(cnt+tam,cnt+1)+Ke(2,3);
            K(cnt+1,cnt+tam+1) = K(cnt+1,cnt+tam+1)+Ke(3,1);
            K(cnt+1,cnt+tam) = K(cnt+1,cnt+tam)+Ke(3,2);
            K(cnt+1,cnt+1) = K(cnt+1,cnt+1)+Ke(3,3);
            
            cnt = cnt+1;
    end
end

%condicoes de contorno para o potencial magnetico (valores conhecidos)
res_b = zeros(tam^2,1);
vpm = ones(tam^2,1);
cnt = 1;
for m = 1:tam
    for n = 1:tam
        if ((m-(tam+1)/2)^2+(n-(tam+1)/2)^2 >= r^2)
            vpm(cnt) = 0;
        elseif(m==(ceil(tam/2)+10/h)&&(n==(ceil(tam/2)-6/h)||(n==(ceil(tam/2)+6/h))))||...
                (m==(ceil(tam/2)+14/h)&&(n==(ceil(tam/2)-4/h)||(n==(ceil(tam/2)+4/h))))
            vpm(cnt) = -mu_ar*imax/(2*pi*rc);
        end
        cnt = cnt + 1;
    end
end

%res_b e a matriz de carregamentos para o vetor potencial magnetico
%mesma ideia do algoritmo anterior
for m = 1:tam^2
    for n = 1:tam^2
        if (vpm(n) ~= 1)&& m == n
            for k = 1:tam^2
                res_b(k) = res_b(k) - K(k,n)*vpm(n);
            end
            K(m,:) = 0;
            K(:,n) = 0;
            K(m,n) = 1;
            res_b(m) = vpm(m);
        end
    end
end

%encontra o potencial magnetico na forma de array (tam^2, 1)
vpm_aux = linsolve(K, res_b);

%passa o potencial magnetico para o dominio 2D
cnt = 1;
Az = zeros(tam);
for m = 1:tam
    for n = 1:tam
        Az(m,n) = vpm_aux(cnt);
        cnt = cnt + 1;
    end
end

%campo eletrico, densidade de fluxo e campo magnetico
%calcula tais variaveis para cada elemento triangular (interpolacao dos
%nohs)
cnt = 1;
Ex = zeros(2*(tam-1));
Ey = Ex;
Bx = Ex;
By = Ex;
Hx = Ex;
Hy = Ex;
for m = 1:tam-1
    for n = 1:tam-1            
            b11 = elementos1(m, n, 4) - elementos1(m, n, 6);
            b21 = elementos1(m, n, 6) - elementos1(m, n, 2);
            b31 = elementos1(m, n, 2) - elementos1(m, n, 4); 
            c11 = elementos1(m, n, 5) - elementos1(m, n, 3); 
            c21 = elementos1(m, n, 1) - elementos1(m, n, 5); 
            c31 = elementos1(m, n, 3) - elementos1(m, n, 1); 
            
            Ae = (b11*c21 - b21*c11)/2;
            
            Ex(2*m-1,2*n-1) = -1/(2*Ae)*(V(m,n)*b11+V(m,n+1)*b21+V(m+1,n)*b31);
            Ey(2*m-1,2*n-1) = -1/(2*Ae)*(V(m,n)*c11+V(m,n+1)*c21+V(m+1,n)*c31);
            By(2*m-1,2*n-1) = 1/(2*Ae)*(Az(m,n)*b11+Az(m,n+1)*b21+Az(m+1,n)*b31);
            Bx(2*m-1,2*n-1) = 1/(2*Ae)*(Az(m,n)*c11+Az(m,n+1)*c21+Az(m+1,n)*c31);
            
            if(m >= ceil(tam/2))
                mu = mu_ar;
            else
                mu = mu_solo;
            end      
            Hx(2*m-1,2*n-1) = Bx(2*m-1,2*n-1)/mu;
            Hy(2*m-1,2*n-1) = By(2*m-1,2*n-1)/mu;
            
            b12 = elementos2(m, n, 4) - elementos2(m, n, 6); 
            b22 = elementos2(m, n, 6) - elementos2(m, n, 2); 
            b32 = elementos2(m, n, 2) - elementos2(m, n, 4); 
            c12 = elementos2(m, n, 5) - elementos2(m, n, 3); 
            c22 = elementos2(m, n, 1) - elementos2(m, n, 5); 
            c32 = elementos2(m, n, 3) - elementos2(m, n, 1); 
            
            Ae = (b12*c22 - b22*c12)/2;
            
            Ex(2*m,2*n) = -1/(2*Ae)*(V(m+1,n+1)*b12+V(m+1,n)*b22+V(m,n+1)*b32);
            Ey(2*m,2*n) = -1/(2*Ae)*(V(m+1,n+1)*c12+V(m+1,n)*c22+V(m,n+1)*c32);
            By(2*m,2*n) = 1/(2*Ae)*(Az(m+1,n+1)*b12+Az(m+1,n)*b22+Az(m,n+1)*b32);
            Bx(2*m,2*n) = 1/(2*Ae)*(Az(m+1,n+1)*c12+Az(m+1,n+1)*c22+Az(m,n+1)*c32);
            
            Hx(2*m,2*n) = Bx(2*m,2*n)/mu;
            Hy(2*m,2*n) = By(2*m,2*n)/mu;
    end
end
            
%aplica condicoes de contorno na fronteira
for n = 1:2:2*(tam-1)
    m = floor(tam/2);
    Ex(m, n+1) = Ex(m+1, n);
    Ey(m, n+1) = sigma_ar*Ey(m+1, n)/sigma_solo;
    By(m, n+1) = By(m+1, n);
    Bx(m, n+1) = mu_ar*Bx(m+1,n)/mu_solo;
    Hx(m, n+1) = Hx(m+1, n);
    Hy(m, n+1) = mu_solo*Hy(m+1,n)/mu_ar;
end

%plota os graficos para V, Az, E, B e H, nessa ordem
close all;
figure(1);
contour(V, 40);
hTitle = title(['Potencial elétrico (kV), h = ', num2str(h)]);
set(hTitle,'FontSize',40);
xlabel('x');
ylabel('y');
set(gca,'FontSize',40);
colormap jet;
colorbar('east');

figure(2);
contour(Az, 40);
hTitle = title(['Potencial vetor magnético (Wb/m), h = ', num2str(h)]);
set(hTitle,'FontSize',40);
xlabel('x');
ylabel('y');
set(gca,'FontSize',40);
colormap jet;
colorbar('east');

figure(3);
x = 1:0.5:tam-0.5;
y = 1:0.5:tam-0.5;
H = quiver(x, y, Ex, Ey, 3);
H.Color = [0.6350 0.0780 0.1840];
hTitle = title(['Campo elétrico (kV/m), h = ', num2str(h)]);
set(hTitle,'FontSize',40);
xlabel('x');
ylabel('y');
set(gca,'FontSize',40);
xlim([1 tam]);
ylim([1 tam-1]);


figure(4);
x = 1:0.5:tam-0.5;
y = 1:0.5:tam-0.5;
quiver(x, y, Bx, -By, 3);
hTitle = title(['Densidade de fluxo magnético (T), h = ', num2str(h)]);
set(hTitle,'FontSize',40);
xlabel('x');
ylabel('y');
set(gca,'FontSize',40);
xlim([1 tam]);
ylim([1 tam-2]);


figure(5);
x = 1:0.5:tam-0.5;
y = 1:0.5:tam-0.5;
H = quiver(x, y, Hx, -Hy, 3);
H.Color = 'black';
hTitle = title(['Campo magnético (A/m), h = ', num2str(h)]);
set(hTitle,'FontSize',40);
xlabel('x');
ylabel('y');
set(gca,'FontSize',40);
xlim([1 tam]);
ylim([1 tam-2]);