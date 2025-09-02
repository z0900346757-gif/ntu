%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                                                         %
%       This code is used for non-parabolic/parabolic band fitting.       %
%                                                                         %
%                                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 注意：QE算出來的倒晶格單位為2pi/alat，要計算的話要轉成1/m，公式為 原值 * 2 * pi / alat / Bohr2Ang / ang2m
% 若是Wannier的話，單位為1/ang，轉成1/m只需要 原值 / ang2m

clear
close all
clc

%% Input / Fitting parameter


% Basic inputs
Material        = 'ZnO';        % Material name
Valley          = 'G';           % Valley name (K: first valley, Q: second valley)
% dEKQ            = 0.2622;        % Unit : eV, Used if "Valley=Q"
FitValence      = 0;             % Set 1 for valence band fitting
iband           = 27;            % band index of conduction/valence band



% Fitting parameters
non_parabolic   = 1;             % Set 1 to fit with non-parabolic equation
alpha           = 0.1;            % Used if non-parabolic == 1, 1.3
%mx              = 0.066601;          % X direction effective mass, 0.13
%my              = 0.06601;           % Y direction effective mass
mx              = 0.226601;          % X direction effective mass, 0.13
my              = 0.226601; 
kmesh           = 200;           % Mesh btw two symmetry points (ex: G-K or K-M)

% Data storing control
store_m_data    = 1;             % Set 1 to store effective mass and alpha data
store_E_data    = 1;             % Set 1 to store 2D Ek data for further SR fitting


% plot control
xmin     = -0.5;
xmax     =  0.5;
ymin     =  -1;
ymax     =  3;
dk       =  0;                   % (修正項) 用於對準 dft data 和 fitting data 的k軸原點，若無必要設為零即可。
plot_2d  =  0;                   % Set 1 to plot 2D band


%% Fundamental/Converting Parameters

e             = 1.6e-19;
me            = 9.11e-31;       % kg
hbar          = 1.05457e-34;     % m^2*kg/s
ang2bohr      = 1.8897259886;  % Å to bohr
bohr2m        = 5.2918e-11;    % bohr to meters
hartree2eV    = 27.2114;       % hartree to eV, 原本的alpha單位是1/hartree
ang2m         = 1e-10;
alat          = 6.122713;     %輸入scf.out的celldm
Ev            = 10.53+0.065;        %輸入band位置

% me, hbar & e = 1 in atomic unit, use this will get hartree energy
e             = 1.6e-19;
me            = 9.11e-31;
hbar          = 1.05457e-34;
mx            = mx*me;
my            = my*me;
pointX        = kmesh*2;
pointY        = kmesh*2;


%% Load file and post-processing DFT raw data

% === 改成讀當前路徑 ===
Filename = [Material,'_',Valley,'.dat'];
Data = importdata(Filename);


% Normalized and plot
tmp1       = find(Data(:,1)==min(Data(:,1)));
tmp2       = find(Data(:,1)==max(Data(:,1)));
tmp3       = [tmp1 tmp2];
tmp4       = tmp3(iband,:);
num_band   = size(tmp3,1);
dft_energy = Data(tmp4(1):tmp4(2),2);   %取出目標band的能量
dft_kpath  = Data(tmp3(1,1):tmp3(1,2),1);  % 取出k點

% Y軸

dft_energy(:) = dft_energy(:) - Ev;


% dft_kpath(:)  = dft_kpath(:)/(ang2m) ; % alat轉成1/m for hse
dft_kpath(:)  = dft_kpath(:) * 2 * pi / alat / bohr2m; % alat轉成1/m for lda pbe

%% Calculate Band Curve
% LDA、PBE
 k1 = 0.9107* 2 * pi / alat / bohr2m; % 1/ang轉1/m  
 k2 = 1.5774* 2 * pi / alat / bohr2m; %
 k3 = 1.8877* 2 * pi / alat / bohr2m; %

%NEW
%k1 = 1.2874* 2 * pi / alat / bohr2m; % 1/ang轉1/m  
%k2 = 1.9540* 2 * pi / alat / bohr2m; %
%k3 = 2.5314* 2 * pi / alat / bohr2m; %
% HSE
% k1 = 1.72596/(ang2m);% 1/ang轉1/m
% k2 = 2.85276/(ang2m);% 1/ang轉1/m
% k3 = 3.45192/(ang2m);% 1/ang轉1/m


num_point = 100;
k_points_left = linspace(k1, k2, num_point);
k_points_right = linspace(k2, k3, num_point);


for i=1:num_point
    if non_parabolic == 1
       % E(1+αE)=ℏ^2(kx^2/(2*mx)+ky^2/(2*my))         
       beta_left         = hbar^2*((k_points_left(i)-k2)^2/(2*mx))/ e;
       EK_2D_left(i) = -1/(2*alpha) + 1/2 * sqrt((1/alpha)^2 + 4*beta_left/alpha);
       beta_right         = hbar^2*((k_points_right(i)-k2)^2/(2*my))/ e;
       EK_2D_right(i) = -1/(2*alpha) + 1/2 * sqrt((1/alpha)^2 + 4*beta_right/alpha);
    else
       % E = ℏ^2(kx^2/(2*mx)+ky^2/(2*my))
       EK_2D_left(i) = hbar^2*((k_points(i)-k2)^2/(2*mx))/ e;  
    end
end

%% Plot band

figure(1)
set(gcf,'Position',[800 100 400 400])

if plot_2d == 1
    
    set(gcf,'Position',[0 0 1400 800])
    subplot(1,2,1);
        
    EK_2D_plot(:,:,1) = kx;
    EK_2D_plot(:,:,2) = ky;
    EK_2D_plot(:,:,3) = EK_2D;   
    
    surf(EK_2D_plot(:,:,2),EK_2D_plot(:,:,1),EK_2D_plot(:,:,3));
    xlabel('k_y ( 2\pi/a )', 'FontSize' ,18);
    ylabel('k_x ( 2\pi/a )', 'FontSize' ,18);
    zlabel('Energy (eV)', 'FontSize' ,16);      
    xlim([-0.7 0.7])
    ylim([-0.7 0.7])
    set(gca,'BoxStyle','full','Box','on')
    ax = gca;
    ax.LineWidth = 2;
    daspect([1  1  8])
    title('2D Band', 'FontSize' ,30);


    subplot(1,2,2);
    grid on
    hold on
    box  on
    ax = gca;
    ax.LineWidth = 2;

    if FitValence == 0
        plot(EK_1D_x(1:pointX/2+1,1),EK_1D_x(1:pointX/2+1,2), '-ro', 'LineWidth',2, 'MarkerEdgeColor','k', 'MarkerFaceColor',[.49 1 .63], 'MarkerSize',3.5);
        plot(EK_1D_y(pointY/2+1:end,1),EK_1D_y(pointY/2+1:end,2), '-bo', 'LineWidth',2, 'MarkerEdgeColor','k', 'MarkerFaceColor',[.49 1 .63], 'MarkerSize',3.5);
    elseif FitValence == 1
        plot(EK_1D_x(1:pointX/2+1,1),-EK_1D_x(1:pointX/2+1,2), '-ro', 'LineWidth',2, 'MarkerEdgeColor','k', 'MarkerFaceColor',[.49 1 .63], 'MarkerSize',3.5);
        plot(EK_1D_y(pointY/2+1:end,1),-EK_1D_y(pointY/2+1:end,2), '-bo', 'LineWidth',2, 'MarkerEdgeColor','k', 'MarkerFaceColor',[.49 1 .63], 'MarkerSize',3.5);
    end
    
    plot(dft_kpath(:)+dk,dft_energy(:), '-go', 'LineWidth',2, 'MarkerEdgeColor','k', 'MarkerFaceColor',[.49 1 .63], 'MarkerSize',3.5);
    daspect([3  2  1])
    
    ylabel('Energy (eV)', 'FontSize' ,16);
    
    if FitValence == 0
        legend('Fitted data X', 'Fitted data Y', 'DFT data','Location','North');
    elseif FitValence == 1
        legend('Fitted data X', 'Fitted data Y', 'DFT data','Location','South');
    end 
    
    xlim([xmin xmax]);
    ylim([ymin ymax]);
    title('Curve Fitting', 'FontSize' ,30);

else
    
    grid on
    hold on
    box  on
    ax = gca;
    ax.LineWidth = 2;

    if FitValence == 0
        plot(k_points_left,EK_2D_left, '-ro', 'LineWidth',2, 'MarkerEdgeColor','k', 'MarkerFaceColor',[.49 1 .63], 'MarkerSize',3.5);
        plot(k_points_right,EK_2D_right, '-bo', 'LineWidth',2, 'MarkerEdgeColor','k', 'MarkerFaceColor',[.49 1 .63], 'MarkerSize',3.5);
    elseif FitValence == 1
        plot(k_points,EK_2D, '-ro', 'LineWidth',2, 'MarkerEdgeColor','k', 'MarkerFaceColor',[.49 1 .63], 'MarkerSize',3.5);
    end    
    
    plot(dft_kpath(:)+dk,dft_energy(:), '-go', 'LineWidth',2, 'MarkerEdgeColor','k', 'MarkerFaceColor',[.49 1 .63], 'MarkerSize',3.5);
    
    ylabel('Energy (eV)', 'FontSize' ,16);
    
    if FitValence == 0
        legend('Fitted data left','Fitted data right',  'DFT data','Location','North');
    elseif FitValence == 1
        legend('Fitted data left','Fitted data right', 'DFT data','Location','South');
    end
    
    ylim([ymin ymax]);
    title('Curve Fitting', 'FontSize' ,30);
    hold off

    
end
