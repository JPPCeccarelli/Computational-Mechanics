function [w2x_d, w2y_d, w2x_e, w2y_e, tensao_Ax, tensao_Ay, tensao_b, D, A, V] = newmark(M, K, C, t_final, dt)

%Dados iniciais
beta = 1/4;
gamma = 1/2;

%matrizes de deslocamento, velocidade, aceleracao e forca
D = zeros(42,t_final/dt+1);
V = zeros(42,t_final/dt+1);
A = zeros(42,t_final/dt+1);
F = zeros(42, 1);

%loop principal
for i = 1:t_final/dt
    F(38,1) = -16000*sin(2*pi*i/t_final);
    F(41,1) = F(38,1);
    if i*dt >= 2 && i*dt <= 8
        F(13,1) = 10000;
        F(19,1) = 10000;
        F(25,1) = 10000;
        F(31,1) = 10000;
        F(37,1) = 10000;
    else
        F(13,1) = 0;
        F(19,1) = 0;
        F(25,1) = 0;
        F(31,1) = 0;
        F(37,1) = 0;
    end   
    
    %Valores atuais para D, V, A (para nao confundir na implementacao das
    %equacoes)
    D_a = D(:,i);
    V_a = V(:,i);
    A_a = A(:,i);
    F1(:,1) = F(:,1)-C*(V_a+dt*(1-gamma)*A_a)-K*(D_a+dt*V_a+(dt^2)*(1-2*beta)*A_a/2);
    A(:,i+1) = ((M+dt*gamma*C+(dt^2)*beta*K)^(-1))*F1(:,1);
    D(:,i+1) = D_a+dt*V_a+(dt^2)*((1-2*beta)*A_a+2*beta*A(:,i+1))/2;
    V(:,i+1) = V_a+dt*((1-gamma)*A_a+gamma*A(:,i+1));
    
    
end

%corrige de acordo com cond de contorno
D(1,:) = 0;
D(2,:) = 0;
D(5,:) = 0;
V(1,:) = 0;
V(2,:) = 0;
V(5,:) = 0;
A(1,:) = 0;
A(2,:) = 0;
A(5,:) = 0;

%deslocamentos em D e E
%cria arrays para registrar o deslocamento nos pontos
w2x_d = zeros(1,t_final/dt);
w2x_e = w2x_d;
w2y_d = w2x_d;
w2y_e = w2x_d;
tensao_Ax = w2x_d;
tensao_Ay = w2x_d;
tensao_b = w2x_d;

%formulas usando polinomios de Hermite
for i = 1:t_final/dt+1
    
    %h = L/2 = 1 nas tensões
    L = 2;
    w2x_d(1, i) = D(16,i)*0.5+D(18,i)*(L-2*((0.5*L)^2)/L+((0.5*L)^3)/L^2)+D(22,i)*0.5+D(24,i)*((-(0.5*L)^2)/L+((0.5*L)^3)/L^2);
    w2y_d(1, i) = D(17,i)*0.5+D(18,i)*(L-2*((0.5*L)^2)/L+((0.5*L)^3)/L^2)+D(23,i)*0.5+D(24,i)*((-(0.5*L)^2)/L+((0.5*L)^3)/L^2);
    tensao_Ax(1, i) = (D(38,i)*(-6/L^2)+D(39,i)*(-4/L)+D(41,i)*6/L^2+D(42,i)*(-2/L))*105*10^9;
    tensao_Ay(1, i) = -(D(37,i)*(-6/L^2)+D(39,i)*(-4/L)+D(31,i)*6/L^2+D(33,i)*(-2/L))*105*10^9;
    
    L = 3/atan(3/0.5);
    w2x_e(1, i) = D(3,i)*(L-2*((0.5*L)^2)/L+((0.5*L)^3)/L^2)+D(7,i)*0.5+D(9,i)*((-(0.5*L)^2)/L+((0.5*L)^3)/L^2);
    w2y_e(1, i) = D(3,i)*(L-2*((0.5*L)^2)/L+((0.5*L)^3)/L^2)+D(8,i)*0.5+D(9,i)*((-(0.5*L)^2)/L+((0.5*L)^3)/L^2);
    
    theta = atan(3/3.5);
    L = 3/cos(theta);
     
    %Obs: h = 0
    tensao_b(1, i) = (-((D(7,i)*cos(theta)+D(8,i)*sin(theta))*(-6/L^2+12*0.5*L/L^3)+...
        D(9,i)*(-4/L+6*0.5*L/L^3)+(D(16,i)*cos(theta)+D(17,i)*sin(theta))*(6/L^2-12*0.5*L/L^3)...
        +D(18,i)*(-2/L+6*0.5*L/L^2))*105*10^9)*L/14;
end

end