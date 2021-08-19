function [f] = f(y)
%f = | f1 = Omega1  |
%    | f2 = Omega2  |
%    | f3 = Omega1' |
%    | f4 = Omega2' |
% Retorna o vetor com as derivadas primeiras correspondentes às 4
% variáveis de estado do vetor y, montado da seguinte maneira:
%y = | y1 = Theta 1  |
%    | y2 = Theta 2  |
%    | y3 = Omega 1  |
%    | y4 = Omega 2  |
theta1 = y(1);
theta2 = y(2);
omega1 = y(3);         
omega2 = y(4);
phi = theta1 - theta2; % Variavel auxiliar utilizada para relacoes trigonometricas, de modo a
                       % simplificar a manipulação algebrica e isolar
                       % Omega1' e Omega2'
                       
%constantes
g = 9.81;
l1 = 2;
l2 = 2.5;
l2e = 1.8;
m1 = 450;
m2 = 650; 
F1 = -0.5*m1*g; 
F2 = -0.5*m2*g;
uIz = 2.7;
R = 0.3;
xd = 80/3.6;
vel = xd;

%equacoes
a0 = (l1^2)*l2*R*(m2*cos(2*phi)-2*m1-m2);
a1 = (l1^2)*l2*R*m2*sin(2*phi);
a2 = 2*l1*(l2^2)*R*m2*sin(phi);
a3 = -2*l2*uIz*vel;
a4 = -2*l1*uIz*vel*cos(phi);
a5 = -R*l1*(l2e*F2*sin(theta1 - 2*theta2)+2*sin(theta1)*(F1*l2+(l2e*F2)/2));

b0 = (l2^2)*R*m2;
b1 = -l1*l2*R*m2*cos(phi);
b2 = l1*l2*R*m2*sin(phi);
b3 = -uIz*vel;
b4 = l2e*sin(theta2)*R*F2;

f3 = (a1*omega1^2+a2*omega2^2+a3*omega1+a4*omega2+a5)/a0;
% Acima, encontra-se a equacao do omega1' isolado
f4 = (b1*f3+b2*(omega1^2)+b3*omega2+b4)/b0;
% Acima, encontra-se a equacao do omega2' isolado
f1 = y(3);
f2 = y(4);
f = [f1 f2 f3 f4];
end