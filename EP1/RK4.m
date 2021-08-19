close all;

h = [1, 0.5, 0.1, 0.01, 0.001]; %Vetor com todos os passos de teste
precision = [0, 1, 1, 2, 3]; %Precisão casas decimais do comando sprintf

for k = 1:length(h) %Resolve as equações para cada passo no vetor h
    y0 = [0 0 0.4 -0.1]; %Vetor com as condições iniciais
    t0 = 0; % Tempo inicial da simulacao
    tf = 60; % Tempo final da simulacao
    N = (tf - t0)/h(k); % Número de passos
    y = y0; %criando o vetor y
    t = t0; %criando vetor t
    p0 = [0.4 -0.1 0 0]; %valores iniciais da plotagem
    p = p0; %cria vetor de plotagem

    for i = 1:N % Loop para executar as iterações do método de Runge Kutta de 4° ordem
        k1 = f(y(i,:));  % K1 = f(Yi)
        k2 = f(y(i,:)+(h(k)/2)*k1); % K2 = f(Yi + (h(k)/2)*k1)
        k3 = f(y(i,:)+(h(k)/2)*k2); % K3 = f(Yi + (h(k)/2)*k2)
        k4 = f(y(i,:)+h(k)*k3); % K4 = f(Yi + h(k)*k3)
        p(i,:) = k1; %passa para o vetor de plotagem
        y(i+1,:) = y(i,:) + (h(k)/6).*(k1+2.*k2+2.*k3+k4); % Y(i+1) = Y(i) + (h(k)/6)/(k1 + 2*k2+ 2*k3 + k4)
        t(i+1) = t(i) + h(k); %Atualizando o tempo de simulacao
    end

    y(end,:) = []; %elimina a ultima linha de y para ter o mesmo tamanho de p
    y(:,end) = []; %elimina a ultima coluna de y para ter o mesmo tamanho de p
    t(end) = []; %elimina a ultima linha de t para ter o mesmo tamanho de p

    %passa os valores de theta1 e theta2 para p
    for i = 1:N-1
        p(i, 5) = y(i, 1);
        p(i, 6) = y(i, 2);
    end

    Plot('RK4', k, p, t, precision, h); %chama a função Plot
end