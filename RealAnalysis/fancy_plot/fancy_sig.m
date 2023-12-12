x = 20:84;
fsize = 16;
y_tick = [6:1:9];
x_tick = [20:10:80,84];

t = tiledlayout(2,2);

nexttile
load("D:\1mat\real_data_analysis\output\online_output\n1000\inference_result\1sig_index_all.mat")
sig_each = sum(sig_index_all(6:70,:)')';
plot(x,sig_each,'-ob','LineWidth',1.5,'MarkerSize',5)
set(gca,'linewidth',1.5,'FontSize',fsize)
axis([20,90,6,9])
set(gca,'YTick',y_tick);
set(gca,'XTick',x_tick);
xlabel('$\textrm{Batch index } m$','Interpreter','latex','FontSize',fsize)
ylabel('$\textrm{Number of significant covariates}$','Interpreter','latex','FontSize',fsize)
title('$n_1=1000$','Interpreter','latex','FontSize',fsize)

clear sig_index_all
clear sig_each

nexttile
load("D:\1mat\real_data_analysis\output\online_output\n3000\inference_result\1sig_index_all.mat")
sig_each = sum(sig_index_all(6:70,:)')';
plot(x,sig_each,'-ob','LineWidth',1.5,'MarkerSize',5)
set(gca,'linewidth',1.5,'FontSize',fsize)
axis([20,90,6,9])
set(gca,'YTick',y_tick);
set(gca,'XTick',x_tick);
xlabel('$\textrm{Batch index } m$','Interpreter','latex','FontSize',fsize)
ylabel('$\textrm{Number of significant covariates}$','Interpreter','latex','FontSize',fsize)
title('$n_1=3000$','Interpreter','latex','FontSize',fsize)


clear sig_index_all
clear sig_each

nexttile
load("D:\1mat\real_data_analysis\output\online_output\n5000\inference_result\1sig_index_all.mat")
sig_each = sum(sig_index_all(6:70,:)')';
plot(x,sig_each,'-ob','LineWidth',1.5,'MarkerSize',5)
set(gca,'linewidth',1.5,'FontSize',fsize)
axis([20,90,6,9])
set(gca,'YTick',y_tick);
set(gca,'XTick',x_tick);
xlabel('$\textrm{Batch index } m$','Interpreter','latex','FontSize',fsize)
ylabel('$\textrm{Number of significant covariates}$','Interpreter','latex','FontSize',fsize)
title('$n_1=5000$','Interpreter','latex','FontSize',fsize)

clear sig_index_all
clear sig_each

nexttile
load("D:\1mat\real_data_analysis\output\online_output\n10000\inference_result\1sig_index_all.mat")
sig_each = sum(sig_index_all(6:70,:)')';
plot(x,sig_each,'-ob','LineWidth',1.5,'MarkerSize',5)
set(gca,'linewidth',1.5,'FontSize',fsize)
axis([20,90,6,9])
set(gca,'YTick',y_tick);
set(gca,'XTick',x_tick);
xlabel('$\textrm{Batch index } m$','Interpreter','latex','FontSize',fsize)
ylabel('$\textrm{Number of significant covariates}$','Interpreter','latex','FontSize',fsize)
title('$n_1=10000$','Interpreter','latex','FontSize',fsize)

t.TileSpacing = 'compact';
%exportgraphics(gcf,filename,'ContentType','vector')