clear; clc; close all;

%% Input
filename = 'ZnO_PBE.dat.gnu';   % bands.x 輸出檔
iband    = 27;                  % <<== ★你想分析的 band 編號 (可以改)
alat     = 3.24;                % Å, QE 的 a 參數
bohr2m   = 5.2918e-11;
ang2m    = 1e-10;

Label = {'\Gamma','M','K','\Gamma','A','L','H','A'};
ymin = -100; ymax = 10;

%% ---- 讀 dat.gnu 檔 ----
fid = fopen(filename);
bands = {};
cur_band = [];
while ~feof(fid)
    tline = fgetl(fid);
    if isempty(tline)
        if ~isempty(cur_band)
            bands{end+1} = cur_band;
            cur_band = [];
        end
    else
        vals = sscanf(tline, '%f %f');
        cur_band = [cur_band; vals'];
    end
end
fclose(fid);
if ~isempty(cur_band)
    bands{end+1} = cur_band;
end

nbnd = length(bands);
k_raw = bands{1}(:,1);
Eband = bands{iband}(:,2);

%% ---- 對齊參考能量線 ----
VBM = max(Eband);  % 對齊到最高點為0
Eband = Eband - VBM;

%% ---- k 軸校正 ----
seg_points = [40,40,40,40,40,40,1];
cum_points = [0, cumsum(seg_points)];
k_ticks = bands{1}(cum_points+1,1);
shift_target = 1.357;
shift_amount = shift_target - k_ticks(4);
k_plot = k_raw + shift_amount;
k_ticks_plot = k_ticks + shift_amount;

%% ---- 自動找 valley（局部極小值）----
valley_indices = [];
for i = 2:length(Eband)-1
    if Eband(i) < Eband(i-1) && Eband(i) < Eband(i+1)
        valley_indices(end+1) = i;
    end
end

% 找出最小能量值
valley_energies = Eband(valley_indices);
min_valley_energy = min(valley_energies);

% 容許些微誤差（浮點數容忍範圍）
tolerance = 1e-3;
true_valley_indices = valley_indices(abs(valley_energies - min_valley_energy) < tolerance);

fprintf('共找到 %d 個 valley:\n', length(true_valley_indices));
for idx = true_valley_indices
    fprintf('  k = %.4f, E = %.4f eV\n', k_plot(idx), Eband(idx));
end

%% ---- 畫圖 ----
figure; hold on; box on;
set(gcf,'Position',[100 100 650 400])
set(gca,'LineWidth',2)

% 畫所有 bands
for ib = 1:nbnd
    plot(k_plot, bands{ib}(:,2) - max(bands{iband}(:,2)), 'Color', [0.8 0.8 0.8]);
end

% 畫你分析的 band
plot(k_plot, Eband, 'b-', 'LineWidth', 2, 'DisplayName', sprintf('Band #%d', iband));

% 標出 valley 點
plot(k_plot(true_valley_indices), Eband(true_valley_indices), 'ro', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Valley');

% 高對稱線
for point = k_ticks_plot
    plot([point point],[ymin ymax],'--k','LineWidth',1,'HandleVisibility','off');
end

% 標籤
xticks(k_ticks_plot);
xticklabels(Label);
set(gca,'XTickLabel',Label,'FontName','Times','FontSize',16)

xlim([min(k_ticks_plot) max(k_ticks_plot)])
ylim([ymin ymax])
ylabel('Energy (eV)','FontWeight','bold','Fontname','Times New Roman','FontSize',18)
title(sprintf('Band #%d Valley Detection', iband), 'FontWeight','bold','Fontname','Times New Roman','FontSize',18)
legend('show','Location','best','FontSize',12)
