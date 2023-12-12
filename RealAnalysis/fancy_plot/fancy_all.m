x = 20:84;
fsize = 16;
y_tick = [0,0.5,1];
x_tick = [20:10:80,84];
lp = 'southwest';
t = tiledlayout(4,2);

load("D:\1mat\real_data_analysis\output\online_output\n1000\inference_result\1p_all_s0.mat")
load("D:\1mat\real_data_analysis\output\online_output\n1000\inference_result\1p_all_s1.mat")
s0_T1 = p_all_s0(6:70,1);
s0_T2 = p_all_s0(6:70,3);
s1_T1 = p_all_s1(6:70,1);
s1_T2 = p_all_s1(6:70,3);

nexttile
plot(x,s0_T1,'-xk','LineWidth',1.5,'MarkerSize',7)
hold on;
plot(x,s0_T2,'--ob','LineWidth',1.5,'MarkerSize',5)
set(gca,'linewidth',1.5,'FontSize',fsize)
axis([20,90,0,1])
set(gca,'YTick',y_tick);
set(gca,'XTick',x_tick);
legend('$T_1$','$T_2$','Interpreter','latex','Location',lp,'Fontsize',fsize,'NumColumns',2)
xlabel('$\textrm{Batch index } m$','Interpreter','latex','FontSize',fsize)
ylabel('$p\textrm{-value}$','Interpreter','latex','FontSize',fsize)
title('($n_1=1000,\sigma_\zeta^2=0$)','Interpreter','latex','FontSize',fsize)

nexttile
plot(x,s1_T1,'-xk','LineWidth',1.5,'MarkerSize',7)
hold on;
plot(x,s1_T2,'--ob','LineWidth',1.5,'MarkerSize',5)
set(gca,'linewidth',1.5,'FontSize',fsize)
axis([20,90,0,1])
set(gca,'YTick',y_tick);
set(gca,'XTick',x_tick);
legend('$T_1$','$T_2$','Interpreter','latex','Location',lp,'Fontsize',fsize,'NumColumns',2)
xlabel('$\textrm{Batch index } m$','Interpreter','latex','FontSize',fsize)
ylabel('$p\textrm{-value}$','Interpreter','latex','FontSize',fsize)
title('($n_1=1000,\sigma_\zeta^2=1$)','Interpreter','latex','FontSize',fsize)

clear p_all_s0
clear p_all_s1
clear s0_T1
clear s0_T2
clear s1_T1
clear s1_T2

load("D:\1mat\real_data_analysis\output\online_output\n3000\inference_result\1p_all_s0.mat")
load("D:\1mat\real_data_analysis\output\online_output\n3000\inference_result\1p_all_s1.mat")
s0_T1 = p_all_s0(6:70,1);
s0_T2 = p_all_s0(6:70,3);
s1_T1 = p_all_s1(6:70,1);
s1_T2 = p_all_s1(6:70,3);

nexttile
plot(x,s0_T1,'-xk','LineWidth',1.5,'MarkerSize',7)
hold on;
plot(x,s0_T2,'--ob','LineWidth',1.5,'MarkerSize',5)
set(gca,'linewidth',1.5,'FontSize',fsize)
axis([20,90,0,1])
set(gca,'YTick',y_tick);
set(gca,'XTick',x_tick);
legend('$T_1$','$T_2$','Interpreter','latex','Location',lp,'Fontsize',fsize,'NumColumns',2)
xlabel('$\textrm{Batch index } m$','Interpreter','latex','FontSize',fsize)
ylabel('$p\textrm{-value}$','Interpreter','latex','FontSize',fsize)
title('($n_1=3000,\sigma_\zeta^2=0$)','Interpreter','latex','FontSize',fsize)

nexttile
plot(x,s1_T1,'-xk','LineWidth',1.5,'MarkerSize',7)
hold on;
plot(x,s1_T2,'--ob','LineWidth',1.5,'MarkerSize',5)
set(gca,'linewidth',1.5,'FontSize',fsize)
axis([20,90,0,1])
set(gca,'YTick',y_tick);
set(gca,'XTick',x_tick);
legend('$T_1$','$T_2$','Interpreter','latex','Location',lp,'Fontsize',fsize,'NumColumns',2)
xlabel('$\textrm{Batch index } m$','Interpreter','latex','FontSize',fsize)
ylabel('$p\textrm{-value}$','Interpreter','latex','FontSize',fsize)
title('($n_1=3000,\sigma_\zeta^2=1$)','Interpreter','latex','FontSize',fsize)

clear p_all_s0
clear p_all_s1
clear s0_T1
clear s0_T2
clear s1_T1
clear s1_T2

load("D:\1mat\real_data_analysis\output\online_output\n5000\inference_result\1p_all_s0.mat")
load("D:\1mat\real_data_analysis\output\online_output\n5000\inference_result\1p_all_s1.mat")
s0_T1 = p_all_s0(6:70,1);
s0_T2 = p_all_s0(6:70,3);
s1_T1 = p_all_s1(6:70,1);
s1_T2 = p_all_s1(6:70,3);

nexttile
plot(x,s0_T1,'-xk','LineWidth',1.5,'MarkerSize',7)
hold on;
plot(x,s0_T2,'--ob','LineWidth',1.5,'MarkerSize',5)
set(gca,'linewidth',1.5,'FontSize',fsize)
axis([20,90,0,1])
set(gca,'YTick',y_tick);
set(gca,'XTick',x_tick);
legend('$T_1$','$T_2$','Interpreter','latex','Location',lp,'Fontsize',fsize,'NumColumns',2)
xlabel('$\textrm{Batch index } m$','Interpreter','latex','FontSize',fsize)
ylabel('$p\textrm{-value}$','Interpreter','latex','FontSize',fsize)
title('($n_1=5000,\sigma_\zeta^2=0$)','Interpreter','latex','FontSize',fsize)

nexttile
plot(x,s1_T1,'-xk','LineWidth',1.5,'MarkerSize',7)
hold on;
plot(x,s1_T2,'--ob','LineWidth',1.5,'MarkerSize',5)
set(gca,'linewidth',1.5,'FontSize',fsize)
axis([20,90,0,1])
set(gca,'YTick',y_tick);
set(gca,'XTick',x_tick);
legend('$T_1$','$T_2$','Interpreter','latex','Location',lp,'Fontsize',fsize,'NumColumns',2)
xlabel('$\textrm{Batch index } m$','Interpreter','latex','FontSize',fsize)
ylabel('$p\textrm{-value}$','Interpreter','latex','FontSize',fsize)
title('($n_1=5000,\sigma_\zeta^2=1$)','Interpreter','latex','FontSize',fsize)

clear p_all_s0
clear p_all_s1
clear s0_T1
clear s0_T2
clear s1_T1
clear s1_T2

load("D:\1mat\real_data_analysis\output\online_output\n10000\inference_result\1p_all_s0.mat")
load("D:\1mat\real_data_analysis\output\online_output\n10000\inference_result\1p_all_s1.mat")
s0_T1 = p_all_s0(6:70,1);
s0_T2 = p_all_s0(6:70,3);
s1_T1 = p_all_s1(6:70,1);
s1_T2 = p_all_s1(6:70,3);

nexttile
plot(x,s0_T1,'-xk','LineWidth',1.5,'MarkerSize',7)
hold on;
plot(x,s0_T2,'--ob','LineWidth',1.5,'MarkerSize',5)
set(gca,'linewidth',1.5,'FontSize',fsize)
axis([20,90,0,1])
set(gca,'YTick',y_tick);
set(gca,'XTick',x_tick);
legend('$T_1$','$T_2$','Interpreter','latex','Location',lp,'Fontsize',fsize,'NumColumns',2)
xlabel('$\textrm{Batch index } m$','Interpreter','latex','FontSize',fsize)
ylabel('$p\textrm{-value}$','Interpreter','latex','FontSize',fsize)
title('($n_1=10000,\sigma_\zeta^2=0$)','Interpreter','latex','FontSize',fsize)

nexttile
plot(x,s1_T1,'-xk','LineWidth',1.5,'MarkerSize',7)
hold on;
plot(x,s1_T2,'--ob','LineWidth',1.5,'MarkerSize',5)
set(gca,'linewidth',1.5,'FontSize',fsize)
axis([20,90,0,1])
set(gca,'YTick',y_tick);
set(gca,'XTick',x_tick);
legend('$T_1$','$T_2$','Interpreter','latex','Location',lp,'Fontsize',fsize,'NumColumns',2)
xlabel('$\textrm{Batch index } m$','Interpreter','latex','FontSize',fsize)
ylabel('$p\textrm{-value}$','Interpreter','latex','FontSize',fsize)
title('($n_1=10000,\sigma_\zeta^2=1$)','Interpreter','latex','FontSize',fsize)

t.TileSpacing = 'compact';
%exportgraphics(gcf,filename,'ContentType','vector')