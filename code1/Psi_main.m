% 
% 
%     l = 30;  % Psi vector length
%     M=10;
%     N=5;
% %     h_mn = sqrt(1/2)*(randn(1,1)+1i*randn(1,1)); % h_mn is plural
% %     h1 = sqrt(1/2)*(randn(1,l)+1i*randn(1,l)); % h1 is a complex vector
% %     h2 = sqrt(1/2)*(randn(1,l)+1i*randn(1,l));  % h2 is a complex vector
%     hmn=sum(h_mn,2).*1000;
%     h1=diag(Theta)';
%     h2=diag(Theta)';
% 
%     % Initialize Psi vector
%     Psi_initial = ones(1, l);
% 
%     % Run gradient descent
%     [Psi_solution, Psi_history, obj_history] = gradient_descent(Psi_initial, hmn, h1, h2);
% 
%     % 绘制图表
%     figure;
%     set(gcf,'position',[0 0 1000 300]);
%     set(gca,'FontSize',9,'Fontname', 'Times New Roman');
% save('Psi_val.mat','Psi_value')
%     subplot(1, 2, 1);
%     set(gcf, 'figsize', [20, 7]);
% 设置子图的宽度和高度
    figure;
%     Psi_history=Psi_history';
    Psi_value=Psi_history;
    colors = ['b','g','r','c','m','y','k'];
    markers = ['o', '+', '*', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h'];
     for i=1:7
         plot(Psi_value(:,i),'LineWidth',1.5,'MarkerSize', 3, 'Marker', markers(mod(i,length(markers))+1));
         hold on;
     end
    legend('group1:1,4,6,21,23,25,27,30', 'group2:2,11,20,28', 'group3:3', 'group4:5,9,13,18,22,26,29','group5:16', 'group6:7,8,12,15,19,24','group7:10,14,17' ,'Location', 'North', 'Orientation', 'horizontal','NumColumns',2);
    xlim([0,50]);ylim([0.8,1.35]);

    
    xlabel('Iteration','Fontname', 'Times New Roman','FontSize',14);
    ylabel('Value of \Psi elements','Fontname', 'Times New Roman','FontSize',14);
    set(gca,'FontSize',12,'FontName', 'Times New Roman'); % 更改坐标轴的字体
    grid on;
%     title('\Psi Over Iterations','Fontname', 'Times New Roman');
    fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print('-painters','-dpdf','-r600','Psi_value_1.pdf')    %前面不管，最后一个参数是pdf名字
%     figure;
% %     subplot(1, 2, 2);
% 
    plot(obj_history,'-o','LineWidth',1.5,'MarkerSize', 4); xlim([0,50]);
    xlabel('Iteration','Fontname', 'Times New Roman');
    ylabel('Objective Function Value','Fontname', 'Times New Roman');
        legend('Algorithm 3,L=30','Location','northeast');
        set(gca,'FontSize',14,'FontName', 'Times New Roman'); % 更改坐标轴的字体
        grid on;
%     title('Objective Function Value Over Iterations','Fontname', 'Times New Roman');
   fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print('-painters','-dpdf','-r600','Obj.pdf')    %前面不管，最后一个参数是pdf名字
