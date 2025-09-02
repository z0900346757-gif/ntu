clear; clc; close all;

%% File names
file_list = {'ZnO_lda.dat.gnu','ZnO_pbe.dat.gnu'};

%% ==== NEW path: Γ → A → H → K → Γ → M → L → H ====
Label = {'\Gamma','A','H','K','\Gamma','M','L','H'};

%% Y-axis range
ymin = -20;
ymax = 10;

%% ---- KPOINTS segment setting (Γ–A–H–K–Γ–M–L–H) ----
seg_points = [40, 40, 40, 40, 40, 40, 40];  % 7 segments
cum_points = [0, cumsum(seg_points)];
k_ticks_index = cum_points + 1;             % indices for Γ,A,H,K,Γ,M,L,H

%% Color & line style
color_list = {'b','r'};       % LDA=blue, PBE=red
style_list_together = {'-','--'};   % LDA=solid, PBE=dashed (together)
style_list_separate = {'-','-'};    % both solid when separate
label_list = {'LDA','PBE'};   % Legend labels

%% Switch: turn=true → same figure, turn=false → separate figures
turn = true;

% Handles for legend
lda_handle = [];
pbe_handle = [];

if turn
    figure; hold on; box on;
    set(gcf,'Position',[100 100 650 400])
    set(gca,'LineWidth',2)
    ylabel('Energy (eV)','FontWeight','bold','Fontname','Times New Roman','FontSize',18)
    title('ZnO Band Structure (LDA vs PBE)','FontWeight','bold','Fontname','Times New Roman','FontSize',18)
end

for f = 1:length(file_list)
    filename = file_list{f};
    
    %% ---- Read file ----
    fid = fopen(filename);
    bands = {};
    cur_band = [];
    while ~feof(fid)
        tline = fgetl(fid);
        if isempty(tline)
            if ~isempty(cur_band)
                bands{end+1} = cur_band; %#ok<SAGROW>
                cur_band = [];
            end
        else
            vals = sscanf(tline, '%f %f');
            cur_band = [cur_band; vals']; %#ok<AGROW>
        end
    end
    fclose(fid);
    if ~isempty(cur_band)
        bands{end+1} = cur_band;
    end

    nbnd = length(bands);
    num_per_band = size(bands{1},1);
    band = zeros(num_per_band, nbnd);
    for ib = 1:nbnd
        band(:,ib) = bands{ib}(:,2);
    end
    
    %% Find VBM
    allE = band(:);
    VBM = max(allE(allE < 0));
    
    %% ---- k axis ----
    k_ticks = bands{1}(k_ticks_index,1);  
    Gamma2_old = k_ticks(5);   % 第二個 Γ 現在是第 5 個節點
    shift_target = 1.357;
    shift_amount = shift_target - Gamma2_old;
    k_plot = bands{1}(:,1) + shift_amount;
    k_ticks_plot = k_ticks + shift_amount;
    
    %% ---- Shift bands ----
    shift_val = 17.246;   % adjust if needed
    band_shifted = band - (VBM + shift_val);
    
    %% ---- Plot ----
    if ~turn
        figure(f); hold on; box on;
        set(gcf,'Position',[200+400*f 200 650 400])
        set(gca,'LineWidth',2)
        ylabel('Energy (eV)','FontWeight','bold','Fontname','Times New Roman','FontSize',18)
        title([label_list{f} ' Band Structure'],'FontWeight','bold','Fontname','Times New Roman','FontSize',18)
    end
    
    % Pick style depending on turn
    if turn
        style_use = style_list_together{f};
    else
        style_use = style_list_separate{f};
    end
    
    % Draw bands, save representative line for legend
    if f == 1
        lda_handle = plot(k_plot, band_shifted(:,1), ...
            'Color',color_list{f}, 'LineStyle',style_use, 'LineWidth',2);
        for ib = 2:nbnd
            plot(k_plot, band_shifted(:,ib), ...
                'Color',color_list{f}, 'LineStyle',style_use, 'LineWidth',2);
        end
    else
        pbe_handle = plot(k_plot, band_shifted(:,1), ...
            'Color',color_list{f}, 'LineStyle',style_use, 'LineWidth',2);
        for ib = 2:nbnd
            plot(k_plot, band_shifted(:,ib), ...
                'Color',color_list{f}, 'LineStyle',style_use, 'LineWidth',2);
        end
    end
    
    % High-symmetry lines (only once in turn=true, always in separate)
    if (turn && f == 1) || (~turn)
        for point = k_ticks_plot
            plot([point point],[ymin ymax],'--k','LineWidth',1,'HandleVisibility','off');
        end
        xticks(k_ticks_plot);
        xticklabels(Label);
        set(gca,'XTickLabel',Label,'FontName','Times','FontSize',16)
        xlim([min(k_ticks_plot) max(k_ticks_plot)])
        ylim([ymin ymax])
        % VBM=0 line
        plot([min(k_ticks_plot) max(k_ticks_plot)],[0 0],'--r','LineWidth',1.2)
    end
    
    %% Show positions in console
    fprintf('\n[%s] High-symmetry k-points:\n', filename);
    for i = 1:length(Label)
        fprintf('%s: x = %.3f\n', Label{i}, k_ticks_plot(i));
    end
end


%% Legend
if turn
    legend([lda_handle pbe_handle], ...
           {'LDA (blue solid)','PBE (red dashed)'}, ...
           'Location','SouthEast','FontSize',12)
else
    if exist('lda_handle','var') && isgraphics(lda_handle)
        legend(lda_handle, {'LDA'}, 'Location','SouthEast','FontSize',12)
    end
    if exist('pbe_handle','var') && isgraphics(pbe_handle)
        legend(pbe_handle, {'PBE'}, 'Location','SouthEast','FontSize',12)
    end
end
