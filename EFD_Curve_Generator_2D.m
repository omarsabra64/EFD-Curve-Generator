clc; clear; close all;

%% --- 1. USER INPUTS (Interactive) ---
disp('--- EFD Parameter Generator ---');
n = input('Enter number of harmonics (n): ');
a0 = input('Enter desired X center (a0): ');
c0 = input('Enter desired Y center (c0): ');

%% --- 2. GENERATE RANDOM SHAPE ---
% We generate a shape centered at (0,0) first to get the shape coefficients
N_points = 1000;
t = linspace(0, 2*pi, N_points);

% Start with a base circle
r = 1.5 * ones(size(t)); 

% Add random "wobbles" to make it a unique shape
% We use low frequencies (2, 3, 4) to keep it smooth and closed
for k_rand = 2:10
    noise_mag = 0.3 * rand(); % Random magnitude
    phase = 2*pi * rand();    % Random phase
    r = r + noise_mag * cos(k_rand * t + phase);
end

% Convert to Cartesian (Centered at 0,0)
x_raw = r .* cos(t);
y_raw = r .* sin(t);

% Force row vectors for calculation
x_raw = reshape(x_raw, 1, []);
y_raw = reshape(y_raw, 1, []);

%% --- 3. CALCULATE COEFFICIENTS (a, b, c, d) ---
t_idx = 1:N_points;
theta = t_idx * (2*pi / N_points);

a = zeros(1, n); 
b = zeros(1, n);
c = zeros(1, n); 
d = zeros(1, n);

for k = 1:n
    % Standard EFD formulas
    a(k) = (2 / N_points) * sum(x_raw .* cos(k * theta));
    b(k) = (2 / N_points) * sum(x_raw .* sin(k * theta));
    c(k) = (2 / N_points) * sum(y_raw .* cos(k * theta));
    d(k) = (2 / N_points) * sum(y_raw .* sin(k * theta));
end

%% --- 4. PLOT (Visual Check) ---
figure('Name', 'Generated EFD Shape', 'Color', 'w');
% Plot the shape shifted to the user's desired a0, c0
plot(x_raw + a0, y_raw + c0, 'b-', 'LineWidth', 2);
hold on;
plot(a0, c0, 'r+', 'MarkerSize', 10, 'LineWidth', 2); % Mark Center
grid on; axis equal;
title(['Generated Shape (n=' num2str(n) ')']);
xlabel('X'); ylabel('Y');
legend('Shape', 'Center (a0, c0)');

%% --- 5. GENERATE COPIABLE OUTPUT ---
fprintf('\n%% --- COPY THE BLOCK BELOW ---\n');
fprintf('a0 = %.2f; c0 = %.2f; n = %d;\n', a0, c0, n);

print_array('a', a);
print_array('b', b);
print_array('c', c);
print_array('d', d);

fprintf('%% ----------------------------\n');

% Helper function
function print_array(var_name, arr)
    fprintf('%s = [', var_name);
    for i = 1:length(arr)
        if i == length(arr)
            fprintf('%.4f', arr(i)); % Last element: no comma
        else
            fprintf('%.4f, ', arr(i)); % Comma separator
        end
    end
    fprintf('];\n');
end