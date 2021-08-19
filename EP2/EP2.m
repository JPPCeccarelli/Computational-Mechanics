%%
delta = 0.375;                                  %escolhe valor para divisoes delta x = delta y
lambda = 1.85;                                  %valor do lambda para sobrerrelaxacao                            
converg = 0.01;                                 %epsilon de convergencia
psi = calcula_psi(delta, lambda, converg);      %funcao que calcula funcao de corrente
[u, v] = velocidades(psi, delta);               %funcao que calcula componentes da velocidade do vento
pressao = calc_pressao(u, v, delta);            %funcao para calcular pressao
[ptx, pty] = pressao_telhado(pressao, delta);   %calcula componentes da pressao ao redor do telhado
forca = calc_forca(delta, pty);                 %calcula forca vertical exercida no telhado

close all
x = linspace(0,36,36/delta + 1);
y = linspace(0,24,24/delta + 1);
[X,Y] = meshgrid(x,y);

%plotagem de figuras
figure;
[C,h] = contour(X,Y,psi,'ShowText','on');
clabel(C,h,'FontSize',20);
hTitle = title("Função corrente \psi do escoamento com \Delta x = " + delta, 'Interpreter', 'tex');
set(hTitle,'FontSize',40);
xlabel('x');
ylabel('y');
set(gca,'FontSize',40);

figure;
mesh(X, Y, psi);
hTitle = title("Função corrente \psi do escoamento com \Delta x = " + delta, 'Interpreter', 'tex');
set(hTitle,'FontSize',40);
xlabel('x');
ylabel('y');
zlabel('\psi');
set(gca,'FontSize',40);

figure;
waterfall(X, Y, psi);
hTitle = title("Função corrente \psi do escoamento com \Delta x = " + delta, 'Interpreter', 'tex');
set(hTitle,'FontSize',40);
xlabel('x');
ylabel('y');
zlabel('\psi');
set(gca,'FontSize',40);

figure;
q = quiver(X,Y,u,v);
q.LineWidth = 2;
q.AutoScaleFactor = 0.7;
hTitle = title("Vetores de velocidade absoluta com \Delta x = " + delta, 'Interpreter', 'tex');
set(hTitle,'FontSize',40);
xlabel('x');
ylabel('y');
set(gca,'FontSize',40);

figure;
[C,h] = contour(X,Y,pressao,'ShowText','on');
clabel(C,h,'FontSize',15);
hTitle = title("Variação de pressão no domínio com \Delta x = " + delta, 'Interpreter', 'tex');
set(hTitle,'FontSize',40);
xlabel('x');
ylabel('y');
set(gca,'FontSize',40);

figure;
mesh(X, Y, pressao);
hTitle = title("Diferença de pressão no domínio com \Delta x = " + delta, 'Interpreter', 'tex');
set(hTitle,'FontSize',40);
xlabel('x');
ylabel('y');
zlabel('Diferença de pressão (Pa)');
set(gca,'FontSize',40);

figure;
q = quiver(ptx,pty);
q.LineWidth = 2;
hTitle = title("Variação de pressão ao redor do telhado com \Delta x = " + delta, 'Interpreter', 'tex');
set(hTitle,'FontSize',40);
set(gca,'FontSize',40);

%%
%Obtencao da distribuicao de temperaturas
delta = 0.375;
lambda = 1.15;
converg = 0.0001;
temp = calcula_temp(delta, lambda, converg, u, v);

x = linspace(0,36,36/delta + 1);
y = linspace(0,24,24/delta + 1);
[X,Y] = meshgrid(x,y);

figure;
[C,h] = contourf(X,Y,temp,'ShowText','on');
clabel(C,h,'FontSize',20)
hTitle = title("Distribuição de Temperatura no domínio com \Delta x = " + delta, 'Interpreter', 'tex');
xlabel('x')
ylabel('y')
set(gca,'FontSize',40);
set(hTitle,'FontSize',40);

%%
delta = 3/8;
fluxo = calcula_fluxo(delta, temp, u, v); %funcao que calcula fluxo de calor saindo do telhado
taxa_calor = calc_taxa_calor(delta,fluxo); %funcao que calcula a taxa de calor total