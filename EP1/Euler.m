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

    for i = 1:N %loop para calcular as iterações no metodo de Euler
        f_ponto = f(y(i,:)); % Derivadas primeiras calculadas na iteração i]
        p(i,:) = f_ponto; %passa para o vetor de plotagem
        y(i+1,:) = y(i,:)+ h(k)*f_ponto; % Y(i+1) = Y(i) + h(k)*f(Y(i))
        t(i+1) = t(i) + h(k); % Atualizando tempo de simulação
    end

    y(end,:) = []; %elimina a ultima linh(k)a de y para ter o mesmo tamanho de p
    y(:,end) = []; %elimina a ultima coluna de y para ter o mesmo tamanho de p
    t(end) = []; %elimina a ultima linh(k)a de t para ter o mesmo tamanho de p

    %passa os valores de theta1 e theta2 para p
    for i = 1:N-1
        p(i, 5) = y(i, 1);
        p(i, 6) = y(i, 2);
    end
    Plot('Euler', k, p, t, precision, h); %chama a função Plot
end