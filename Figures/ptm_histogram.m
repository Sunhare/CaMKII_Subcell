%% Histogram
cd '/Users/satolab/Documents/Projects/CaMKII_Subcell/'

orig_data = load('ave.txt');
% ptm_actual_data = load('ave_sptm_1000.txt');
% ptm_data = load('ave_ptm_SATO.txt');
ptm_data = load('ave_ptm.txt');

nryr = orig_data(:,1);
vp = orig_data(:,2);
Jmax = orig_data(:,3);
ci = orig_data(:,4);
cs = orig_data(:,5);
cp = orig_data(:,6);


nryr_ptm = ptm_data(:,1);
vp_ptm = ptm_data(:,2);
Jmax_ptm = ptm_data(:,3);
ci_ptm = ptm_data(:,4);
cs_ptm = ptm_data(:,5);
cp_ptm = ptm_data(:,6);


%Plot #1
figure('Position', [0 0 1920/3 1080/2])
set(gcf, 'Renderer', 'Painters')
set(gcf, 'color','w')
h = histogram(nryr,12,'FaceColor', [1 0 0], 'EdgeColor', 'black', 'FaceAlpha',1, 'LineWidth', 1)
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


%Plot #3
figure('Position', [1920*2/3 0 1920/3 1080/2])
set(gcf, 'Renderer', 'Painters')
set(gcf, 'color','w')
h2 = histogram(cp_ptm*1000,18,'FaceColor', [1 0 0], 'EdgeColor', 'black', 'FaceAlpha',1, 'LineWidth', 1);
% set(h, 'barwidth', 0.9)
xlabel('[Ca]_C_l_e_f_t (nM)')
ylabel('# of CRUs')
set(gca, 'FontName',  'arial', 'FontSize', 50)
set(gca,'linewidth',1)

box off
ylim([0 200])
xlim([0 1000])

%Plot #4
figure('Position', [0 1920*2/3 1920/3 1080/2])
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


%% CJSR plot
ptm_data = load('/Users/satolab/Documents/Projects/CaMKII_Subcell/results/ave_sptm.txt');
nryr_ptm = ptm_data(:,2);
cjsr_ptm = ptm_data(:,7);

ec50SR = 0.45;
MaxSR = 15; MinSR = 1;

kCaSR = MaxSR - (MaxSR-MinSR)./(1+(ec50SR./(cjsr_ptm.*1e-3)).^2.5);
koSRCa = 10./kCaSR;

figure
h = histogram(cjsr_ptm.*1e-3,12,'FaceColor', [1 0 0], 'EdgeColor', 'black', 'FaceAlpha',1, 'LineWidth', 1);
% title("C_J_S_R Distribution")
xlabel("[Ca]_J_S_R (mM)", 'FontSize', 25)
set(gca, 'FontSize', 25)
xlim([0 2])
figure
h = histogram(koSRCa,12,'FaceColor', [1 0 0], 'EdgeColor', 'black', 'FaceAlpha',1, 'FaceColor', 'blue', 'LineWidth', 1);
xlabel("koSRCa", 'FontSize', 25 )
xlim([3 4])
set(gca, 'FontSize', 25)
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

%% 
ptm_data = load('ave_sptm.txt');
figure
nryr_ptm = ptm_data(:,1);
cp_ptm = ptm_data(:,6);
plot(nryr_ptm, cp_ptm,'+k')
title("PTM")

figure
plot(nryr,cp, '+k')
title("Control")

%%
%TODO  Avg [Ca](cleft space) vs cluster size

scatter(nryr, cp,100,'k','+') %Original
% scatter(nryr_ptm, cp_ptm,100,'k','+') %PTM
