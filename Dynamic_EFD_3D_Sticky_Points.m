clc; clear; close all;

%% --- 1. SETUP ---
disp('--- Real-Time 3D EFD Animator ---');
n = input('Enter number of harmonics (n): ');
num_sticky = input('Enter number of sticky points (X): '); 

% Default centers
a0 = 0; c0 = 0; e0 = 0; 

% Animation settings
N_points = 500; 
t = linspace(0, 2*pi, N_points);

% --- STICKY POINTS SETUP ---
% Define the fixed t values for the sticky points (evenly spaced)
t_sticky = linspace(0, 2*pi, num_sticky+1); 
t_sticky(end) = []; % Remove last point to avoid overlap with start
% Generate distinct colors for each point
sticky_colors = hsv(num_sticky); 

% --- NOISE PARAMETERS ---
n_noise = 5;
% Pool size is 10 to ensure randperm doesn't error
xy_mag = 0.1 + 0.2*rand(1, n_noise);    
xy_freq = 2 + randperm(10, n_noise);     
xy_phase = 2*pi*rand(1, n_noise);        

z_mag = 0.2 + 0.3*rand(1, n_noise);      
z_freq = 1 + randperm(10, n_noise);      
z_phase = 2*pi*rand(1, n_noise);         

% Evolution speed
drift_speed = 0.005; 

%% --- 2. FIGURE SETUP ---
f = figure('Name', 'Real-Time EFD Generation', 'Color', 'w');

% Plot Objects
h_plot = plot3(nan, nan, nan, 'b-', 'LineWidth', 1.5); hold on;

% ERROR FIX: Initialize scatter with vectors of NaNs matching num_sticky size
h_sticky = scatter3(nan(1, num_sticky), nan(1, num_sticky), nan(1, num_sticky), ...
                    60, sticky_colors, 'filled'); 

h_center = plot3(a0, c0, e0, 'r+', 'MarkerSize', 10, 'LineWidth', 3);

grid on; axis equal; view(3);
xlabel('X'); ylabel('Y'); zlabel('Z');
xlim([-2.5 2.5]); ylim([-2.5 2.5]); zlim([-1.5 1.5]);

% Stop Button
stop_btn = uicontrol('Style', 'pushbutton', 'String', 'Stop & Copy',...
    'Position', [20 20 100 40], 'Callback', 'set(gcf,''UserData'',1)');

% Set initial state
f.UserData = 0; 

%% --- 3. REAL-TIME LOOP ---
disp('Animating... Press the "Stop & Copy" button on the figure to end.');

% Pre-allocate coefficient arrays
a = zeros(1, n); b = zeros(1, n);
c = zeros(1, n); d = zeros(1, n);
e_coef = zeros(1, n); f_coef = zeros(1, n);

while isvalid(f) && f.UserData == 0
    tic; % Start Timer
    
    % -- A. EVOLVE SHAPE --
    xy_phase = xy_phase + drift_speed;
    z_phase = z_phase + drift_speed;
    
    % Generate Raw Target Points
    r = 1.5 * ones(size(t));
    z_raw = zeros(size(t));
    
    for k = 1:n_noise
        r = r + xy_mag(k) * cos(xy_freq(k) * t + xy_phase(k));
        z_raw = z_raw + z_mag(k) * sin(z_freq(k) * t + z_phase(k));
    end
    x_raw = r .* cos(t);
    y_raw = r .* sin(t);
    
    % -- B. CALCULATE COEFFICIENTS --
    for k = 1:n
        cos_k = cos(k * t);
        sin_k = sin(k * t);
        
        a(k) = (2 / N_points) * sum(x_raw .* cos_k);
        b(k) = (2 / N_points) * sum(x_raw .* sin_k);
        c(k) = (2 / N_points) * sum(y_raw .* cos_k);
        d(k) = (2 / N_points) * sum(y_raw .* sin_k);
        e_coef(k) = (2 / N_points) * sum(z_raw .* cos_k);
        f_coef(k) = (2 / N_points) * sum(z_raw .* sin_k);
    end
    
    % -- C. RECONSTRUCT MAIN CURVE --
    x_rec = a0 * ones(size(t));
    y_rec = c0 * ones(size(t));
    z_rec = e0 * ones(size(t));
    
    for k = 1:n
        cos_k = cos(k * t);
        sin_k = sin(k * t);
        x_rec = x_rec + a(k)*cos_k + b(k)*sin_k;
        y_rec = y_rec + c(k)*cos_k + d(k)*sin_k;
        z_rec = z_rec + e_coef(k)*cos_k + f_coef(k)*sin_k;
    end
    
    % -- D. CALCULATE STICKY POINTS --
    % We compute positions specifically for t_sticky using current coeffs
    x_s = a0 * ones(size(t_sticky));
    y_s = c0 * ones(size(t_sticky));
    z_s = e0 * ones(size(t_sticky));
    
    for k = 1:n
        cos_k_s = cos(k * t_sticky);
        sin_k_s = sin(k * t_sticky);
        x_s = x_s + a(k)*cos_k_s + b(k)*sin_k_s;
        y_s = y_s + c(k)*cos_k_s + d(k)*sin_k_s;
        z_s = z_s + e_coef(k)*cos_k_s + f_coef(k)*sin_k_s;
    end
    
    % -- E. UPDATE PLOTS --
    set(h_plot, 'XData', x_rec, 'YData', y_rec, 'ZData', z_rec);
    set(h_sticky, 'XData', x_s, 'YData', y_s, 'ZData', z_s); % Update dots
    
    % Measure performance
    dt = toc;
    fps = 1/dt;
    title({['Real-Time EFD (n=' num2str(n) ')']; ...
           ['FPS: ' num2str(fps, '%.1f') ' | Calc Time: ' num2str(dt*1000, '%.1f') ' ms']});
    
    drawnow limitrate; 
end

%% --- 4. OUTPUT FINAL STATE ---
if isvalid(f)
    clc;
    fprintf('\n%% --- FINAL CAPTURED STATE ---\n');
    fprintf('a0 = %.2f; c0 = %.2f; e0 = %.2f; n = %d;\n', a0, c0, e0, n);
    print_array('a', a);
    print_array('b', b);
    print_array('c', c);
    print_array('d', d);
    print_array('e', e_coef);
    print_array('f', f_coef);
    fprintf('%% ----------------------------\n');
else
    disp('Figure closed manually. No data copied.');
end

function print_array(var_name, arr)
    fprintf('%s = [', var_name);
    for i = 1:length(arr)
        if i == length(arr)
            fprintf('%.4f', arr(i)); 
        else
            fprintf('%.4f, ', arr(i)); 
        end
    end
    fprintf('];\n');
end