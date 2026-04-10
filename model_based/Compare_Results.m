clear all
close all
clc
%%
Farrage_loadsway = readmatrix('Farrage_loadSway.xlsx');
Farrage_theta4   = readmatrix('Farrage_position.xlsx');
Farrage_dtheta_4 = readmatrix('Farrage_velocity.xlsx');

w_Farrage = [Farrage_loadsway, -(Farrage_theta4(:,2) - Farrage_theta4(1,2)), Farrage_dtheta_4(:,2)];

Iskandar_loadsway = readmatrix('Iskander_loadSway.xlsx');
Iskandar_theta4   = readmatrix('Iskander_position.xlsx');
Iskandar_dtheta_4 = readmatrix('Iskander_velocity.xlsx');

w_Iskandar = [Iskandar_loadsway(1:end-2,:), -(Iskandar_theta4(1:end-2,2) - Iskandar_theta4(1,2)), Iskandar_dtheta_4(1:end-2,2)];

% mssing the target by
% missing target in model based
deg2rad(45 - w_Farrage(end,4))
% missing target in model based
deg2rad(45 - w_Iskandar(end,4))
% percentage of the error
(deg2rad(45 - w_Farrage(end,4))/deg2rad(45 - w_Iskandar(end,4)))*100


%% LaTeX settings
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaulttextinterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

% Create figure with background color
fig = figure('Color', [232/255, 232/255, 232/255]);


% Colors
color_ref   = [0.9290, 0.6940, 0.1250];  % orange
color_opt   = [0, 0.3176, 0.6196];       % blue
color_meas  = [1, 0, 0];                 % red
color_dd    = [0.8500, 0.3250, 0.0980];  % optional for w_dd
color_xline = [0.5, 0.5, 0.5];           % grey
color_Farr  = [0.0, 0.38, 0.18];    % green (similar intensity to the blue)


% Subplot positions (4x1 vertical layout)
nSub = 4;
bottomMargin = 0.08;
topMargin    = 0.05;
verticalSpacing = 0.05;
subplotHeight = (1 - topMargin - bottomMargin - (nSub-1)*verticalSpacing)/nSub;
subplotLeft = 0.12;
subplotWidth = 0.8;

positions = zeros(nSub,4);
for i = 1:nSub
    positions(i,:) = [subplotLeft, bottomMargin + (nSub-i)*(subplotHeight + verticalSpacing), subplotWidth, subplotHeight];
end

% Plot properties
fs_label  = 9;
fs_tick   = 9;
fs_title  = 10;
lw_line   = 1;
ms_marker = 4;

% --- Loop over 4 subplots ---
for k = 1:4
    ax = axes('Position', positions(k,:));
    hold on;

    % Farrage
    plot(w_Farrage(:,1), deg2rad(w_Farrage(:,k+1)), '--', 'LineWidth', lw_line, 'Color', color_opt);
    % measured
    plot(w_Iskandar(:,1), deg2rad(w_Iskandar(:,k+1)), '-', 'LineWidth', lw_line, 'Color', color_meas);
    hold off;

    % Labels
    if k < 4
        set(gca,'XTick', []);
    else
        xlabel('Time [s]', 'FontSize', fs_label, 'Interpreter','latex');
    end
    if k < 3
        ylabel(sprintf('$\\theta_%d$ [rad]', k), 'FontSize', fs_label, 'Interpreter','latex');
    elseif k==3
        ylabel(sprintf('$\\theta_%d$ [rad]', k+1), 'FontSize', fs_label, 'Interpreter','latex');
    end
       
    set(gca,'FontSize', fs_tick, 'XLim', [15,75]);

    if k==1 || k==2
        ylim([-0.03,0.03])
    elseif k==3
        ylim([-0.2,1])
    else
        ylim([-0.02,0.18])
    end
    box on;
    if k==4
        hold on
        h_Farr    = plot(nan, nan, '--', 'Color', color_opt, 'LineWidth', lw_line);
        h_meas  = plot(nan, nan, '-', 'Color', color_meas, 'LineWidth', lw_line);
        hold off
    
    legend([h_Farr, h_meas], ...
           {'Model-based', ...
            'Data-Driven'}, ...
           'Location', 'northoutside', ...
           'Orientation', 'horizontal', ...
           'Interpreter', 'latex', ...
           'FontSize', fs_label, ...
           'Box', 'off');    
        ylabel('$\dot{\theta}_4$ [rad/s]', 'FontSize', fs_label, 'Interpreter','latex');
    end
end

% Adjust figure dimensions
width = 9;  % cm (long and narrow)
height = 8; % cm total
set(fig,'Units','centimeters')
set(fig, 'Position', [0, 0, width, height]);
set(fig, 'PaperSize', [width, height]);
set(fig,'PaperPositionMode','auto');


%% Simulation model validation
experimental_loadsway = readmatrix('Experimental_LoadSway.xlsx');
simulation_loadsway   = readmatrix('Simulation_LoadSway.xlsx');
validation_dtheta_4   = readmatrix('Horizontal_Velocity.xlsx');

% Set LaTeX interpreters
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaulttextinterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

% Colors
color_meas  = [0.0000, 0.3176, 0.6196]; % grey
color_opt   = [1, 0, 0]; % blue
color_vline = [0.5, 0.5, 0.5];         % grey for vertical lines

% Figure
fig = figure('Color',[1 1 1]);

% Paper-format dimensions (cm)
fig_width  = 9.2;
fig_height = 7;
set(fig,'Units','centimeters','Position',[0,0,fig_width,fig_height], ...
    'PaperUnits','centimeters','PaperSize',[fig_width,fig_height]);

% Adjusted subplot layout (2x2) for better spacing
left   = 0.14;  % left margin for left column
bottom = 0.14;  % bottom margin for bottom row
w      = 0.33;  % width of each subplot
h      = 0.30;  % height of each subplot
hgap   = 0.12;  % horizontal gap between subplots
vgap   = 0.12;  % vertical gap between subplots

positions = [
    left,           bottom + h + vgap + 0.02,      w, h;   % top-left (move slightly up)
    left + w + hgap + 0.03, bottom + h + vgap + 0.02, w, h;   % top-right (move right and up)
    left,           bottom,                        w, h;   % bottom-left
    left + w + hgap + 0.03, bottom,                w, h;   % bottom-right (move right)
];

fs_label = 9; % font size for labels

% Loop over the four subplots
for k = 1:4
    ax = subplot('Position', positions(k,:));
    hold on
    
    % Plot simulated (nonparametric) trajectory
    plot(experimental_loadsway(:,1), deg2rad(experimental_loadsway(:,k+1)), '-', 'LineWidth',1,'Color',color_opt);
    
    % Plot measured trajectory
    plot(simulation_loadsway(:,1), deg2rad(simulation_loadsway(:,k+1)), '--', 'LineWidth',1,'Color',color_meas);
    
    % Scatter sampled measured points for visibility
    scatter(simulation_loadsway(1:25:end,1), deg2rad(simulation_loadsway(1:25:end,k+1)), 10, ...
        'MarkerFaceColor', color_meas, 'MarkerEdgeColor', color_meas, 'MarkerFaceAlpha',0.2);
   
    % Axis labels
    if k > 2
        xlabel('Time [s]','FontSize',fs_label,'Interpreter','latex');
    end

    switch k
        case 1, ylabel('$\theta_1$ [rad]','FontSize',fs_label,'Interpreter','latex');
        case 2, ylabel('$\dot{\theta}_1$ [rad/s]','FontSize',fs_label,'Interpreter','latex');
        case 3, ylabel('$\theta_2$ [rad]','FontSize',fs_label,'Interpreter','latex');
        case 4, ylabel('$\dot{\theta}_2$ [rad/s]','FontSize',fs_label,'Interpreter','latex');
    end
    
    % Axis limits
    switch k
        case {1,2}, ylim([-0.005,0.005]);
        case 3,    ylim([-0.005,0.005]);
        case 4,    ylim([-0.005,0.005]);
    end
    xlim([0,20])
    if k == 4
        % Legend on the last subplot
        h_meas = plot(nan,nan,'--','Color',color_meas,'LineWidth',1);
        h_opt  = plot(nan,nan,'-','Color',color_opt,'LineWidth',1);
        h_input  = plot(nan,nan,'-','Color',color_vline,'LineWidth',1);
        
        legend([h_opt , h_meas, h_input], {'Experimental','Simulation', 'Input'}, ...
            'Location','northoutside','Orientation','horizontal', ...
            'Interpreter','latex','FontSize',fs_label,'Box','off');
    end
    set(gca,'FontSize',fs_label)
    box on
    hold off
end

% Create new figure with same height as subplots
fig2 = figure('Color',[1 1 1]);

% Match the same paper size (width identical to the 2-column layout)
fig_width  = 9.2;
fig_height = 3.5;   % same height as one subplot
set(fig2, 'Units','centimeters', 'Position',[0,0,fig_width,fig_height], ...
    'PaperUnits','centimeters','PaperSize',[fig_width,fig_height]);

% Convert normalized margins from previous figure to cm
left_margin_cm = 0.14 * fig_width;  % left margin in cm
right_margin_cm = fig_width - (left_margin_cm + 7.455);  % to get exact 74.55 mm axis width

% Convert margins back to normalized units for 'Position'
left_norm = left_margin_cm / fig_width;
width_norm = 7.455 / fig_width;  % exact width = 74.55 mm
bottom_norm = 0.2;
height_norm = 2.12 / fig_height;

% Create axes with the exact width
ax_wide = axes('Position', [left_norm, bottom_norm, width_norm, height_norm]);
hold on

% Plot
plot(validation_dtheta_4(:,1), deg2rad(validation_dtheta_4(:,2)), ...
    'Color', color_vline, 'LineWidth', 1.2);

xlabel('Time [s]', 'FontSize', fs_label, 'Interpreter','latex');
ylabel('$\dot{\theta}_4$ [rad/s]', 'FontSize', fs_label, 'Interpreter','latex');
ylim([-0.05,0.25])
set(gca, 'FontSize', fs_label);
box on
hold off




