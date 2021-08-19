function Plot(metodo, k, p, t, precision, h)

    fontSize = 36; %Tamanho da fonte do título
    fontSizeAxis = 24; %Tamanho da fonte dos eixos
    lineWidth = 3; %Espessura da linha do gráfico
    
    %Inicia o subplot
    
    figure(k*2-1);
    subplot(3,2,1);
    plot(t,p(:,5),'LineWidth', lineWidth);
    grid on;
    xlabel('$t[s]$', 'interpreter', 'latex','FontSize', fontSize);
    ylabel('$\theta_1[rad]$', 'interpreter', 'latex','FontSize', fontSize);
    set(gca,'FontSize',fontSizeAxis);

    subplot(3,2,2);
    plot(t,p(:,6),'LineWidth', lineWidth);
    grid on;
    xlabel('$t[s]$', 'interpreter', 'latex','FontSize', fontSize);
    ylabel('$\theta_2[rad]$', 'interpreter', 'latex','FontSize', fontSize);
    set(gca,'FontSize',fontSizeAxis);

    subplot(3,2,3);
    plot(t,p(:,1),'LineWidth', lineWidth);
    grid on;
    xlabel('$t[s]$', 'interpreter', 'latex','FontSize', fontSize);
    ylabel('$\dot\theta_1[rad/s]$', 'interpreter', 'latex','FontSize', fontSize);
    set(gca,'FontSize',fontSizeAxis);

    subplot(3,2,4);
    plot(t,p(:,2),'LineWidth', lineWidth);
    grid on;
    xlabel('$t[s]$', 'interpreter', 'latex','FontSize', fontSize);
    ylabel('$\dot\theta_2[rad/s]$', 'interpreter', 'latex','FontSize', fontSize);
    set(gca,'FontSize',fontSizeAxis);

    subplot(3,2,5);
    plot(t,p(:,3),'LineWidth', lineWidth);
    grid on;
    xlabel('$t[s]$', 'interpreter', 'latex','FontSize', fontSize);
    ylabel('$\ddot\theta_1[rad/s^2]$', 'interpreter', 'latex','FontSize', fontSize);
    set(gca,'FontSize',fontSizeAxis);

    subplot(3,2,6);
    plot(t,p(:,4),'LineWidth', lineWidth);
    grid on;
    xlabel('$t[s]$', 'interpreter', 'latex','FontSize', fontSize);
    ylabel('$\ddot\theta_2[rad/s^2]$', 'interpreter', 'latex','FontSize', fontSize);
    set(gca,'FontSize',fontSizeAxis);
    
    sgtitle(sprintf('Simulação para método de %s para passo h = %1.*f', metodo, precision(k), h(k)),'FontSize', fontSize);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

    saveas(gcf,sprintf('C:\\Users\\henri\\Documents\\EP1\\EP1\\ex1\\%sh%1.*fseparado.jpg',metodo, precision(k),h(k))); %Salva o arquivo em .JPG
    
    %Plota todas as curvas no mesmo gráfico
    figure(k*2)
    plot(t, p(:,5), t, p(:,6), t, p(:,1), t, p(:,2), t, p(:,3), t, p(:,4),'LineWidth', lineWidth);%plotagem do grafico
    title(sprintf('Simulação para método de %s para passo h = %1.*f', metodo, precision(k), h(k)),'FontSize', fontSize) %Titulo do grafico
    legend('$\theta_1[rad]$', '$\theta_2[rad]$', '$\dot\theta_1[rad/s]$', '$\dot\theta_2[rad/s]$', '$\ddot\theta_1[rad/s^2]$', '$\ddot\theta_2[rad/s^2]$', 'interpreter', 'latex', 'FontSize', fontSize) %Legenda das curvas
    xlabel('$t[s]$', 'interpreter', 'latex','FontSize', fontSize) % Titulo do eixo x
    ylabel('$[rad]$ | $[rad/s]$ | $[rad/s^2]$','interpreter', 'latex','FontSize', fontSize) % Titulo do eixo y
    set(gca,'FontSize',fontSizeAxis);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(gcf,sprintf('C:\\Users\\henri\\Documents\\EP1\\EP1\\ex1\\%sh%1.*f.jpg', metodo, precision(k),h(k)));  %Salva o arquivo em .JPG
end

