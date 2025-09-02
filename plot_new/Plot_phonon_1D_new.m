% Points of high symmetry
points = {'\Gamma', 'K', 'M', '\Gamma', 'A'};

points_coor = zeros(length(points),1);

mode = 1; % mode = 0 for cm-1, mode = 1 for meV 

% === Load data ===
freq_data = load('./GaAs_band/ZnO3.freq.gp'); % 改成你的路徑
x = freq_data(:,1);

count = 1;
for i = 1:length(x)
    if mod(i,40) == 1
        points_coor(count) = x(i);
        count = count + 1;
    end
end

% === Y data (frequencies) ===
if mode == 0
    y = freq_data(:,2:end);
    ymax = 800;
    y_unit = 'cm^{-1}';
else
    y = freq_data(:,2:end) * 1.239841E-1; % cm^-1 → meV
    ymax = 100;
    y_unit = 'meV';
end

xmin = points_coor(1);
xmax = points_coor(end);
ymin = 0;

% === Plot ===
figure;
hold on;
for i = 1:size(y,2)
    plot(x, y(:,i), 'blue', 'LineWidth', 2);
end

% Lines from points of high symmetry
for i = 1:length(points)
    plot([points_coor(i) points_coor(i)], [ymin ymax], 'Color', [0.5 0.5 0.5]);
end

% Axis properties
xlim([xmin xmax]);
ylim([ymin ymax]);
ax = gca;
box on
set(gcf,'Position',[1 48 400 350])
set(gca,'LineWidth',2)

if mode == 0
    ylabel(['Frequency (' y_unit ')'],'FontWeight','bold','Fontname','Times New Roman','FontSize',24);
    yticks(ymin:100:ymax);
else
    ylabel(['Energy (' y_unit ')'],'FontWeight','bold','Fontname','Times New Roman','FontSize',24);
    yticks(ymin:10:ymax);
end

xticks(points_coor);
xticklabels(points);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',20)

% === Title ===
title('ZnO phonon dispersion','FontWeight','bold','Fontname','Times New Roman','FontSize',24);

% === 計算並顯示總 phonon 數 ===
n_modes = size(y,2);
disp(['Total number of phonon branches = ', num2str(n_modes)]);
text(xmax*0.7, ymax*0.9, ['Total modes = ' num2str(n_modes)], ...
     'FontSize',16,'FontName','Times','FontWeight','bold');

% === 找出最高點 ===
[max_val, idx] = max(y(:));        % 找最大值
[row, col] = ind2sub(size(y), idx); % 對應座標
x_max = x(row);

% 在 console 顯示
disp(['Highest phonon frequency = ' num2str(max_val) ' ' y_unit ...
      ' at x = ' num2str(x_max)]);

% 在圖上標出來
plot(x_max, max_val, 'ro', 'MarkerSize', 8, 'LineWidth', 2); % 紅點
text(x_max, max_val+5, ['Max = ' num2str(max_val, '%.2f') ' ' y_unit], ...
     'FontSize',14,'FontName','Times','FontWeight','bold','Color','red');
