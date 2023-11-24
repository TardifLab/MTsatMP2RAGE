% Plots of data used in qMT_LinearModelVals
% This needs to be separated since not all values are reported in each
% study.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 5 parameters, 4 free ones as R1 is always X. Assume Ra = R1obs = Rb
figure; tiledlayout(2,2);
nexttile;

%% R1 vs M0b
R1 = [ 0.99, 1.72, 0.95, 1.7, 1/0.806, 1/0.873, 1/0.839,1/1.070, 1/1.458,...
    1/1.043, 1/1.064, 1/1.1121, 1/1.088, 1/1.522, 1/1.372, 1/1.407,...
    0.97, 0.71, 1/0.794, 1.04, 0.68, 0.89];
M0b = [0.056, 0.152, 0.072, 0.161,0.133, 0.122, 0.134, 0.155, 0.094,...
    0.153, 0.149, 0.139, 0.1489, 0.0643, 0.0658, 0.088,...
    0.1536, 0.088,  0.113, 0.114, 0.075, 0.125];


x1_line = linspace( 0.5, 1.75, 50);
y1 = polyfit( R1, M0b, 1); y2 = polyval(y1, x1_line);
caption1 = [ 'M_{0B} ≈', num2str(y1(1),  '%.3g'),'·R_{1obs} + ', num2str(y1(2),  '%.3g') ] ;

plot(x1_line, y2,'LineWidth',3,'Color','k'); hold on;
scatter(R1, M0b, 30, [0, 0.4470, 0.7410],'filled')
ax = gca; ax.FontSize = 20; xlim([0.5, 1.75])
xlabel( 'R_{1obs}', 'FontSize', 18, 'FontWeight', 'bold')
ylabel('M_{0B}', 'FontSize', 18, 'FontWeight', 'bold')
text(0.55, 0.162, caption1,"FontSize",16)
ylim([0.05, 0.17]);
nexttile;

%% R1 vs kf
R1 = [ 0.99, 1.72, 0.95, 1.7, 1/0.806, 1/0.873, 1/0.839,...
    1/1.043, 1/1.064, 1/1.1121, 1/1.088, 1/1.522, 1/1.372, 1/1.407,...
    0.97, 0.71, 1/0.794, 1.04, 0.68];
kf = [2.2, 4.6, 2.4, 4.3, 3.3, 3.1, 3.1,...
    2.758, 2.57, 2.196, 1.956, 1.322, 1.072, 1.358, 2.7, 1.57, 5, 1.254, 1.125];

x1_line = linspace( 0.5, 1.75, 50);
y1 = polyfit( R1, kf, 1); y2 = polyval(y1, x1_line);
caption1 = [ 'k_{f} ≈', num2str(y1(1),  '%.3g'),'·R_{1obs} ', num2str(y1(2),  '%.3g') ] ;

plot(x1_line, y2,'LineWidth',3,'Color','k'); hold on;
scatter(R1, kf, 30, [0, 0.4470, 0.7410],'filled')
ax = gca; ax.FontSize = 20; xlim([0.5, 1.75])
xlabel( 'R_{1obs}', 'FontSize', 18, 'FontWeight', 'bold')
ylabel('k_{f}', 'FontSize', 18, 'FontWeight', 'bold')
ylim([0, 5.2])
text(0.55, 4.2, caption1,"FontSize",16)
nexttile;

%% R1 vs T2b
R1 = [ 0.99, 1.72, 0.95, 1.7, 0.97, 0.71];
T2b = [9.7, 11.8, 11.1, 12.3, 9.84, 10.21]*10^(-6);

x1_line = linspace( 0.5, 1.75, 50);
y1 = polyfit( R1, T2b, 1); y2 = polyval(y1, x1_line);
caption1 = [ 'T_{2B} ≈', num2str(y1(1),  '%.3g'),'·R_{1obs} + ', num2str(y1(2),  '%.3g') ] ;

plot(x1_line, y2,'LineWidth',3,'Color','k'); hold on;
scatter(R1, T2b, 30, [0, 0.4470, 0.7410],'filled')
ax = gca; ax.FontSize = 20; xlim([0.5, 1.75])
xlabel( 'R_{1obs}', 'FontSize', 18, 'FontWeight', 'bold')
ylabel('T_{2B}', 'FontSize', 18, 'FontWeight', 'bold')
text(0.55, 12.5e-6, caption1,"FontSize",16)
ylim([9e-6, 13e-6]);
nexttile;

%% R1 vs T2a
R1 = [ 0.99, 1.72, 0.95, 1.7, 1/0.794];
T2a = [55, 31, 56, 37, 54]*10^(-3);

x1_line = linspace( 0.5, 1.75, 50);
y1 = polyfit( R1, T2a, 1); y2 = polyval(y1, x1_line);
caption1 = [ 'T_{2A} ≈', num2str(y1(1),  '%.3g'),'·R_{1obs} + ', num2str(y1(2),  '%.3g') ] ;

plot(x1_line, y2,'LineWidth',3,'Color','k'); hold on;
scatter(R1, T2a, 30, [0, 0.4470, 0.7410],'filled')
ax = gca; ax.FontSize = 20; xlim([0.5, 1.75])
xlabel( 'R_{1obs}', 'FontSize', 18, 'FontWeight', 'bold')
ylabel('T_{2A}', 'FontSize', 18, 'FontWeight', 'bold')
text(0.55, 35e-3, caption1,"FontSize",16)
ylim([0.03,0.075])



set(gcf,'position',[10,400,1100,800])  








