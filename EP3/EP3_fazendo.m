%PMR3401 - EP3
%Alunos: João Pedro Ceccarelli e Henrique Yda Yamamoto

%Dados iniciais
L1 = 2;
L2 = 3;
L3 = 4;
d1i = 0.072;
d1e = 0.08;
d2i = 0.09;
d2e = 0.1;
E = 210*10^9;
ro = 7650;
delta_t = 0.5;
t_final = 20;

%matriz de massa local e global
%dimensão 42 pois são 14 nós que contém 3 variáveis cada
Me = zeros(6);
M = zeros(42);

%matriz de rigidez local e global
Ke = zeros(6);
K = zeros(42);

%matriz de rotacao
R = zeros(6);

phi = [];

%insere os angulos de cada barra da torre em phi
for i = 1:27
    if i==1||i==6||i==12||i==16||i==20||i==24||i==27
        phi(i) = 0;
    elseif i==11||i==13||i==15||i==17||i==19||i==21||i==23||i==25
        phi(i) = pi/2;
    elseif i==14||i==22
        phi(i) = pi/4;
    elseif i==18||i==26
        phi(i) = 3*pi/4;
    elseif i==2||i==7
        phi(i) = atan(3/0.5);
    elseif i==3||i==8
        phi(i) = atan(-3/0.5)+pi;
    elseif i==4
        phi(i) = atan(3/3.5);
    elseif i==5
        phi(i) = atan(-3/3.5)+pi;
    elseif i==10
        phi(i) = atan(-3/2.5)+pi;
    else %9
        phi(i) = atan(3/2.5);
    end
end

%vetor de nós nas vigas, que identifica as extremidades das vigas
nv = [1 2; 1 3; 2 4; 1 4; 2 3; 3 4; 3 5; 4 6; 3 6; 4 5; 5 7; ...
    5 6; 6 8; 5 8; 7 9; 7 8; 8 10; 8 9; 9 11; 9 10; 10 12; 9 12; 11 13; ...
    11 12; 12 14; 12 13; 13 14];

for i = 1:27
    %atualiza dados para cada barra diferente
    if i~=4 && i~=5 && i~=9 && i~=10 && i~=14 && i~=18 && i~=22 && i~=26
        A = pi*(d2e^2 - d2i^2)/4;
        I = pi*(d2e^4 - d2i^4)/64;
        if i==1
            L = L3;
        elseif i==2||i==3||i==7||i==8
            L = L2/sin(phi(i));
        elseif i==6
            L = L2;
        else
            L = L1;
        end
    else
        A = pi*(d1e^2 - d1i^2)/4;
        I = pi*(d1e^4 - d1i^4)/64;
        if i==4||i==5||i==9||i==10
            L = L2/sin(phi(i));
        else
            L = L1/sin(phi(i));
        end
    end
    
    %atualiza matriz de rotacao
    R = [cos(phi(i)) sin(phi(i)) 0 0 0 0;
        -sin(phi(i)) cos(phi(i)) 0 0 0 0;
        0 0 1 0 0 0;
        0 0 0 cos(phi(i)) sin(phi(i)) 0;
        0 0 0 -sin(phi(i)) cos(phi(i)) 0;
        0 0 0 0 0 1];
    
    %Matriz de massa retirada de Timoshenko
    Me = [140 0 0 70 0 0;
        0 156 22*L 0 54 -13*L;
        0 22*L 4*L^2 0 13*L -3*L^2;
        70 0 0 140 0 0;
        0 54 13*L 0 156 -22*L;
        0 -13*L -3*L^2 0 -22*L 4*L^2];
    
    Me = Me*(ro*A/420);
    
    %rotação da matriz local
    Me = R'*Me*R;
    
    %coloca na matriz de massa global: cria uma matriz auxiliar que
    %transfere os valores das matrizes locais para a matriz global
    aux = zeros(42);
    aux(3*(nv(i,1)-1)+1:3*nv(i,1), 3*(nv(i,1)-1)+1:3*nv(i,1)) = Me(1:3, 1:3);
    aux(3*(nv(i,1)-1)+1:3*nv(i,1), 3*(nv(i,2)-1)+1:3*nv(i,2)) = Me(1:3, 4:6);
    aux(3*(nv(i,2)-1)+1:3*nv(i,2), 3*(nv(i,1)-1)+1:3*nv(i,1)) = Me(4:6, 1:3);
    aux(3*(nv(i,2)-1)+1:3*nv(i,2), 3*(nv(i,2)-1)+1:3*nv(i,2)) = Me(4:6, 4:6);
    
    M = M + aux;
    
    %Agora fazendo o mesmo para a matriz de rigidez
    Ke = [E*A/L 0 0 -E*A/L 0 0;
        0 12*E*I/L^3 6*E*I/L^2 0 -12*E*I/L^3 6*E*I/L^2;
        0 6*E*I/L^2 4*E*I/L 0 -6*E*I/L^2 2*E*I/L;
        -E*A/L 0 0 E*A/L 0 0;
        0 -12*E*I/L^3 -6*E*I/L^2 0 12*E*I/L^3 -6*E*I/L^2;
        0 6*E*I/L^2 2*E*I/L 0 -6*E*I/L^2 4*E*I/L];
    
    Ke = R'*Ke*R;
    
    aux = zeros(42);
    aux(3*(nv(i,1)-1)+1:3*nv(i,1), 3*(nv(i,1)-1)+1:3*nv(i,1)) = Ke(1:3, 1:3);
    aux(3*(nv(i,1)-1)+1:3*nv(i,1), 3*(nv(i,2)-1)+1:3*nv(i,2)) = Ke(1:3, 4:6);
    aux(3*(nv(i,2)-1)+1:3*nv(i,2), 3*(nv(i,1)-1)+1:3*nv(i,1)) = Ke(4:6, 1:3);
    aux(3*(nv(i,2)-1)+1:3*nv(i,2), 3*(nv(i,2)-1)+1:3*nv(i,2)) = Ke(4:6, 4:6);
    
    K = K + aux;
end

%funcao que calcula as frequencias de vibracao
[y_x, y_y, y_m, lambda] = modos_de_vibracao(K, M);

%funcao que calcula os deslocamentos
[desc, desb, descf, desbf] = resp_frequencia(y_x, y_y, y_m);

%matriz de amortecimento
C = 0.3*M + 0.03*K;

[wxd, wyd, wxe, wye, t_ax, t_ay, t_b, D, A, V] = newmark(M, K, C, t_final, delta_t);

%plotagem dos graficos
%modos de vibração
close all;
figure;
plot(y_x(:,1),'red','lineWidth',2); hold on;
plot(y_x(:,2),'blue','lineWidth',2); hold on;
plot(y_x(:,3),'magenta','lineWidth',2); hold on;
plot(y_x(:,4),'black','lineWidth',2); hold on;
plot(y_x(:,5),'yellow','lineWidth',2); hold on;
plot(y_x(:,6),'green','lineWidth',2);
xlabel('Nó');
ylabel('Deformação horizontal (m)');
legend('1','2','3','4','5','6');
title('Modos de vibração');

figure;
plot(y_y(:,1),'red','lineWidth',2); hold on;
plot(y_y(:,2),'blue','lineWidth',2); hold on;
plot(y_y(:,3),'magenta','lineWidth',2); hold on;
plot(y_y(:,4),'black','lineWidth',2); hold on;
plot(y_y(:,5),'yellow','lineWidth',2); hold on;
plot(y_y(:,6),'green','lineWidth',2);
xlabel('Nó');
ylabel('Deformação vertical (m)');
legend('1','2','3','4','5','6');
title('Modos de vibração');

figure;
plot(y_m(:,1),'red','lineWidth',2); hold on;
plot(y_m(:,2),'blue','lineWidth',2); hold on;
plot(y_m(:,3),'magenta','lineWidth',2); hold on;
plot(y_m(:,4),'black','lineWidth',2); hold on;
plot(y_m(:,5),'yellow','lineWidth',2); hold on;
plot(y_m(:,6),'green','lineWidth',2);
xlabel('Nó');
ylabel('Deformação angular (rad)');
legend('1','2','3','4','5','6');
title('Modos de vibração');

%deslocamentos dos pontos B e C
figure;
bar(lambda(1:3,:), desc(1:3, :), 0.01);
xlabel('Frequência (Hz)');
ylabel('||u|| (m)');
title('Resposta em frequência - Hermite (Ponto C)');

figure;
bar(lambda(1:3,:), desb(1:3, :), 0.01);
xlabel('Frequência (Hz)');
ylabel('||u|| (m)');
title('Resposta em frequência - Hermite (Ponto B)');

figure;
bar(lambda(1:3,:), descf(1:3, :), 0.01);
xlabel('Frequência (Hz)');
ylabel('||u|| (m)');
title('Resposta em frequência - Funções de forma simples (Ponto C)');

figure;
bar(lambda(1:3,:), desbf(1:3, :), 0.01);
xlabel('Frequência (Hz)');
ylabel('||u|| (m)');
title('Resposta em frequência - Funções de forma simples (Ponto B)');

%deslocamentos
figure;
plot(0:delta_t:t_final, wxd(1,:),'red','lineWidth',2); hold on;
plot(0:delta_t:t_final, wyd(1,:),'blue','lineWidth',2);
xlabel('Tempo (s)');
ylabel('Deslocamento (m)');
legend('Direção X', 'Direção Y');
title('Deslocamentos no ponto D');

figure;
plot(0:delta_t:t_final, wxe(1,:),'green','lineWidth',2); hold on;
plot(0:delta_t:t_final, wye(1,:),'magenta','lineWidth',2);
xlabel('Tempo (s)');
ylabel('Deslocamento (m)');
legend('Direção X', 'Direção Y');
title('Deslocamentos no ponto E');

%tensões
figure;
plot(0:delta_t:t_final, t_ax(1,:),'blue','lineWidth',0.5); hold on;
plot(0:delta_t:t_final, t_ay(1,:),'red','lineWidth',0.5);
xlabel('Tempo (s)');
ylabel('Tensão (Pa)');
legend('Ponto A - direção X', 'Ponto A - direção Y');
title('Tensão');

figure;
plot(0:delta_t:t_final, t_b(1,:),'green','lineWidth',0.5);
xlabel('Tempo (s)');
ylabel('Tensão (Pa)');
legend('Ponto B');
title('Tensão');

%Respostas transientes
figure;
plot(0:delta_t:t_final, D(19,:), 'blue', 'lineWidth',0.5); hold on;
plot(0:delta_t:t_final, V(19,:), 'red', 'lineWidth',0.5); hold on;
plot(0:delta_t:t_final, A(19,:), 'green', 'lineWidth',0.5);
xlabel('Tempo (s)');
ylabel('');
legend('Deslocamento (m)', 'Velocidade (m/s)', 'Aceleração (m/s^2)');
title('Nó 7');

figure;
plot(0:delta_t:t_final, D(20,:), 'blue', 'lineWidth',0.5); hold on;
plot(0:delta_t:t_final, V(20,:), 'red', 'lineWidth',0.5); hold on;
plot(0:delta_t:t_final, A(20,:), 'green', 'lineWidth',0.5);
xlabel('Tempo (s)');
ylabel('');
legend('Deslocamento (m)', 'Velocidade (m/s)', 'Aceleração (m/s^2)');
title('Nó 7 (y)');

figure;
plot(0:delta_t:t_final, D(21,:), 'blue', 'lineWidth',0.5); hold on;
plot(0:delta_t:t_final, V(21,:), 'red', 'lineWidth',0.5); hold on;
plot(0:delta_t:t_final, A(21,:), 'green', 'lineWidth',0.5);
xlabel('Tempo (s)');
ylabel('');
legend('Deslocamento (m)', 'Velocidade (m/s)', 'Aceleração (m/s^2)');
title('Nó 7 (\theta)');

figure;
plot(0:delta_t:t_final, D(7,:), 'blue', 'lineWidth',0.5); hold on;
plot(0:delta_t:t_final, V(7,:), 'red', 'lineWidth',0.5); hold on;
plot(0:delta_t:t_final, A(7,:), 'green', 'lineWidth',0.5);
xlabel('Tempo (s)');
ylabel('');
legend('Deslocamento (m)', 'Velocidade (m/s)', 'Aceleração (m/s^2)');
title('Nó 3');

figure;
plot(0:delta_t:t_final, D(8,:), 'blue', 'lineWidth',0.5); hold on;
plot(0:delta_t:t_final, V(8,:), 'red', 'lineWidth',0.5); hold on;
plot(0:delta_t:t_final, A(8,:), 'green', 'lineWidth',0.5);
xlabel('Tempo (s)');
ylabel('');
legend('Deslocamento (m)', 'Velocidade (m/s)', 'Aceleração (m/s^2)');
title('Nó 3 (y)');

figure;
plot(0:delta_t:t_final, D(9,:), 'blue', 'lineWidth',0.5); hold on;
plot(0:delta_t:t_final, V(9,:), 'red', 'lineWidth',0.5); hold on;
plot(0:delta_t:t_final, A(9,:), 'green', 'lineWidth',0.5);
xlabel('Tempo (s)');
ylabel('');
legend('Deslocamento (m)', 'Velocidade (m/s)', 'Aceleração (m/s^2)');
title('Nó 3 (\theta)');