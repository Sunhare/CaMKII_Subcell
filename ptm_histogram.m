%% Histogram
cd '/Users/satolab/Documents/Projects/CaMKII_Subcell/results'

ptm_data = load('ave.txt');

nryr = ptm_data(:,1);
vp = ptm_data(:,2);
Jmax = ptm_data(:,3);
ci = ptm_data(:,4);
cs = ptm_data(:,5);
cp = ptm_data(:,6);


%Plot #1
figure('Position', [0 0 1920/3 1080/2])
set(gcf, 'Renderer', 'Painters')
set(gcf, 'color','w')
h = histogram(nryr,12,'FaceColor', [1 0 0], 'EdgeColor', 'black', 'FaceAlpha',1, 'LineWidth', 1)
% set(h, 'barwidth', 0.9)
xlabel('# RyRs')
ylabel('# of CRUs')
set(gca, 'FontName',  'arial', 'FontSize', 50)
set(gca,'linewidth',1)

box off
ylim([0 200])

%Plot #2
figure('Position', [1920/3 0 1920/3 1080/2])
set(gcf, 'Renderer', 'Painters')
set(gcf, 'color','w')
h2 = histogram(cp*1000,12,'FaceColor', [1 0 0], 'EdgeColor', 'black', 'FaceAlpha',1, 'LineWidth', 1);
% set(h, 'barwidth', 0.9)
xlabel('[Ca]_C_l_e_f_t (nM)')
ylabel('# of CRUs')
set(gca, 'FontName',  'arial', 'FontSize', 50)
set(gca,'linewidth',1)

box off
ylim([0 200])
xlim([0 1000])

%% I-V Curve
figure
set(gcf, 'Color', 'white')
data_iv = load("0d_results_vc_20.txt");
            
t_iv = data_iv(:,1);
v_iv = data_iv(:,2);
jca_iv = data_iv(:,3);
ica_ave = data_iv(:,4);
jnaca_iv = data_iv(:,6);
x1=-20;
x2 = 480;

subplot 311
plot(t_iv, v_iv, 'LineWidth', 5)
ylabel("Voltage (mV)")
legend("Voltage")
ylim([-90 30])
xlim([x1 x2])
set(gca, 'FontSize', 25)
subplot 312
plot(t_iv, ica_ave, 'b','LineWidth', 5)
legend("I_C_a (ave)")
xlim([x1 x2])
set(gca, 'FontSize', 25)
subplot 313
plot(t_iv, jnaca_iv, 'Color', [0.3, 0.8, 0.5], 'LineWidth',5)
legend("J_N_a_C_a", 'FontSize', 25)
xlim([x1 x2])
set(gca, 'FontSize', 25)