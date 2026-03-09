% --- Parameters ---
m = 0.6;   % slope of regression line
b = 0.5;   % intercept

% Observed point above the line (further away for clarity)
x_obs = 2.5;
y_obs = 3.8;

% OLS projection (vertical)
y_hat_ols = m*x_obs + b;

% TLS projection (perpendicular)
x_hat_tls = (x_obs + m*(y_obs - b)) / (1 + m^2);
y_hat_tls = m*x_hat_tls + b;

% --- Figure setup ---
fig_width_cm  = 12; % width in cm
fig_height_cm = 10;  % height in cm
figure;
set(gcf, 'Units', 'centimeters', 'Position', [2 2 fig_width_cm fig_height_cm]);

hold on;
%axis equal;

% Wider limits for more annotation space
xlim([1.0 5.0]);   % expanded horizontally
ylim([0.5 4.8]);   % expanded vertically

% Regression line
x_line = linspace(0,5,200);
y_line = m*x_line + b;
plot(x_line, y_line, 'k-', 'LineWidth', 1.5);

% Points
plot(x_obs, y_obs, 'ko', 'MarkerFaceColor', 'k');
plot(x_obs, y_hat_ols, 'ks', 'MarkerFaceColor', [0.7 0.7 0.7]);
plot(x_hat_tls, y_hat_tls, 'kd', 'MarkerFaceColor', [0.7 0.7 0.7]);

% OLS vertical distance
plot([x_obs x_obs], [y_obs y_hat_ols], 'r-', 'LineWidth', 2);
mid_ols_y = (y_obs + y_hat_ols)/2;
x_offset = 0.15;
text(x_obs - x_offset, mid_ols_y, {'Ordinary', 'least squares', 'distance'}, ...
    'Color', 'r', 'FontSize', 12, 'Interpreter', 'latex', ...
    'HorizontalAlignment', 'right');

% TLS perpendicular distance
plot([x_obs x_hat_tls], [y_obs y_hat_tls], 'b-', 'LineWidth', 2);
mid_tls_x = (x_obs + x_hat_tls)/2;
mid_tls_y = (y_obs + y_hat_tls)/2;
text(mid_tls_x + 0.25, mid_tls_y + 0.1, {'Orthogonal', 'least squares', 'distance'}, ...
    'Color', 'b', 'FontSize', 12, 'Interpreter', 'latex');

% Right angle marker for TLS
ra_size = 0.1;
dx = x_hat_tls - x_obs;
dy = y_hat_tls - y_obs;
len = sqrt(dx^2 + dy^2);
ux = dx/len; uy = dy/len; % unit vector along TLS
vx = -uy;    vy = ux;     % perpendicular unit vector
p1 = [x_hat_tls, y_hat_tls];
corner = p1 + ra_size * [-ux, -uy];
corner2 = corner + ra_size * [vx, vy];
corner3 = p1 + ra_size * [vx, vy];
plot([corner(1) corner2(1)], [corner(2) corner2(2)], 'k-', 'LineWidth', 1);
plot([corner2(1) corner3(1)], [corner2(2) corner3(2)], 'k-', 'LineWidth', 1);

% Labels for points
text(x_obs + 0.15, y_obs + 0.15, '$\left(w_i, y_i\right)$', ...
    'Interpreter', 'latex', 'FontSize', 16);
text(x_hat_tls + 0.15, y_hat_tls - 0.15, '$\left(x_i, y(\theta,x_i)\right)$', ...
    'Interpreter', 'latex', 'FontSize', 16);

% Axis labels
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 18);

box on;
set(gca, 'FontSize', 16);

% Save as PNG at 300 dpi
%print(gcf, 'OLSillustration.png', '-dpng', '-r300');

hold off;
