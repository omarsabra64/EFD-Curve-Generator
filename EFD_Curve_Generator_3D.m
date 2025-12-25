clc; clear; close all;

%% --- 1. USER INPUTS (Interactive) ---
disp('--- 3D EFD Parameter Generator & Benchmark ---');
n = input('Enter number of harmonics (n): ');
a0 = input('Enter desired X center (a0): ');
c0 = input('Enter desired Y center (c0): ');
e0 = input('Enter desired Z center (e0): ');

%% --- 2. GENERATE RANDOM SHAPE ---
% Generating raw points to fit
N_points = 1000;
t = linspace(0, 2*pi, N_points);

% -- Generate X and Y (Wobbly Circle) --
r = 1.5 * ones(size(t)); 
for k_rand = 2:6
    r = r + 0.3 * rand() * cos(k_rand * t + 2*pi*rand());
end
x_raw = r .* cos(t);
y_raw = r .* sin(t);

% -- Generate Z (Vertical Undulation) --
z_raw = zeros(size(t));
for k_z = 1:5
    z_raw = z_raw + 0.5 * rand() * sin(k_z * t + 2*pi*rand());
end

% Reshape for calculation
x_raw = reshape(x_raw, 1, []);
y_raw = reshape(y_raw, 1, []);
z_raw = reshape(z_raw, 1, []);

%% --- 3. CALCULATE COEFFICIENTS (BENCHMARK 1) ---
fprintf('\nCalculating Coefficients...\n');
tic; % <--- START TIMER 1 (Fitting Speed)

t_idx = 1:N_points;
theta = t_idx * (2*pi / N_points);
a = zeros(1, n); b = zeros(1, n);
c = zeros(1, n); d = zeros(1, n);
e_coef = zeros(1, n); f = zeros(1, n); % Renamed 'e' to 'e_coef' to avoid conflict with exp(1)

for k = 1:n
    % Standard EFD formulas (Discrete Fourier Transform)
    a(k) = (2 / N_points) * sum(x_raw .* cos(k * theta));
    b(k) = (2 / N_points) * sum(x_raw .* sin(k * theta));
    c(k) = (2 / N_points) * sum(y_raw .* cos(k * theta));
    d(k) = (2 / N_points) * sum(y_raw .* sin(k * theta));
    e_coef(k) = (2 / N_points) * sum(z_raw .* cos(k * theta));
    f(k) = (2 / N_points) * sum(z_raw .* sin(k * theta));
end

fit_time = toc; % <--- STOP TIMER 1

%% --- 4. RECONSTRUCTION LOOP (BENCHMARK 2) ---
% Measure how fast we can generate points FROM the coefficients
% (This simulates the computational load on your robot controller)
fprintf('Benchmarking Point Generation...\n');

x_rec = zeros(1, N_points);
y_rec = zeros(1, N_points);
z_rec = zeros(1, N_points);

tic; % <--- START TIMER 2 (Generation Speed)

for i = 1:N_points
    curr_t = theta(i);
    
    % Reset coordinate values
    val_x = a0; val_y = c0; val_z = e0;
    
    % Sum the harmonics
    for k = 1:n
        cos_val = cos(k*curr_t);
        sin_val = sin(k*curr_t);
        
        val_x = val_x + a(k)*cos_val + b(k)*sin_val;
        val_y = val_y + c(k)*cos_val + d(k)*sin_val;
        val_z = val_z + e_coef(k)*cos_val + f(k)*sin_val;
    end
    x_rec(i) = val_x;
    y_rec(i) = val_y;
    z_rec(i) = val_z;
end

gen_time = toc; % <--- STOP TIMER 2

%% --- 5. REPORT PERFORMANCE ---
fprintf('\n--- PERFORMANCE METRICS ---\n');
fprintf('1. Model Fitting Time (Raw -> Coeffs): %.6f seconds\n', fit_time);
fprintf('2. Curve Generation Time (Coeffs -> Points): %.6f seconds\n', gen_time);
fprintf('   - Average time per point: %.8f seconds\n', gen_time/N_points);
fprintf('   - Max Control Frequency:  %.2f kHz\n', (N_points/gen_time)/1000);

%% --- 6. PLOT ---
figure('Name', '3D EFD Shape & Reconstruction', 'Color', 'w');
plot3(x_rec, y_rec, z_rec, 'b-', 'LineWidth', 2); hold on;
plot3(a0, c0, e0, 'r+', 'MarkerSize', 10, 'LineWidth', 3);
grid on; axis equal; view(3);
title(['Generated Shape (n=' num2str(n) ')']);
xlabel('X'); ylabel('Y'); zlabel('Z');
legend('Reconstructed Path', 'Center');

%% --- 7. OUTPUT ---
fprintf('\n%% --- COPY THE BLOCK BELOW ---\n');
fprintf('a0 = %.2f; c0 = %.2f; e0 = %.2f; n = %d;\n', a0, c0, e0, n);
print_array('a', a);
print_array('b', b);
print_array('c', c);
print_array('d', d);
print_array('e', e_coef);
print_array('f', f);
fprintf('%% ----------------------------\n');

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