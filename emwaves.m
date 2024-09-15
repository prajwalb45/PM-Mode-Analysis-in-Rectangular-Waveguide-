clc;
close all;

% Waveguide dimensions
a = 2.286;  % Length in cm in x-direction
b = a/2;    % Length in cm in y-direction
f = 45*10.^9;   % Frequency of operation 45GHz
c = 3*10.^8;    % Velocity of light 

% User input for mode values
m = input('Enter mode value m:');
n = input('Enter mode value n:');    

Amn = 1;    % Particular mode Constant

% Wave propagation in Z-Direction
%********************************%
fc = c*100/2*sqrt((m/a)^2+(n/b)^2);    % Cutoff frequency calculation in GHz
lambda = c*100/fc;              % Wavelength in cm
epsilon = 8.8540e-12;           % Permittivity constant
epsilon_r = 1;                  % Relative Permittivity constant
mu1 = 4*pi*10e-7;               % Permeability constant
mu1_r = 1;                      % Relative Permeability constant
omega = 2*pi*f;                 % Frequency of operation in rad/s
M = 40;                         % Number of points to be plotted

beta = omega*(sqrt(mu1*epsilon));  % Propagation constant
Bx = m*pi/a;    % Beta(x)
By = n*pi/b;    % Beta(y)
Bc = sqrt(Bx^2 + By^2); % Beta(c), cutoff wavenumber
Bz = sqrt(beta^2 - Bc^2);

% Check for mode existence and cutoff frequency
if m == 0 || n == 0
    fprintf(['TM_', num2str(m), num2str(n), ' mode does not exist']);
elseif fc > f
    fprintf(['TM_', num2str(m), num2str(n), ' mode cutoff frequency exceeds frequency of operation; hence mode does not propagate\n']);
    sprintf('The frequency of operation is up to: %0.5g', f)
    sprintf('The cutoff frequency is: %0.5g', fc)
else
    sprintf('The frequency of operation is up to: %0.5g', f)
    sprintf('The cutoff frequency is: %0.5g', fc)

    % Field Pattern plot for Rectangular wave guide for TMmn mode
    % Front View
    x = linspace(0, a, M);
    y = linspace(0, b, M);
    [x, y] = meshgrid(x, y);

    % Placeholder field expressions (replace with actual expressions)
    Ex = Amn * cos(Bx * x) .* sin(By * y) .* exp(-1i * Bz * 0);
    Ey = Amn * sin(Bx * x) .* cos(By * y) .* exp(-1i * Bz * 0);
    Ez = Amn * sin(Bx * x) .* sin(By * y) .* exp(-1i * Bz * 0);
    Hx = Amn * sin(m * pi * x / a) .* cos(n * pi * y / b) .* exp(-1i * Bz * 0);
    Hy = Amn * cos(m * pi * x / a) .* sin(n * pi * y / b) .* exp(-1i * Bz * 0);
    Hz = zeros(size(Hy));  % Assuming Hz is zero for simplicity

    % Plot of TMmn E-Field view
    figure();
    quiver(x, y, real(Ex), real(Ey));
    title(['Plot of front view for TM_', num2str(m), '_', num2str(n), ' E-Field']);
    legend('E-Field');
    xlabel('x-dimension 0 to a');
    ylabel('y-dimension 0 to b=a/2');

    % Plot of TMmn H-Field view
    figure();
    quiver(x, y, real(Hx), real(Hy));
    title(['Plot of front view for TM_', num2str(m), '_', num2str(n), ' H-Field']);
    legend('H-Field');
    xlabel('x-dimension 0 to a');
    ylabel('y-dimension 0 to b=a/2');

    % Plot of TMmn E-Field and H-Field view
    figure();
    quiver(x, y, real(Ex), real(Ey));
    hold on
    quiver(x, y, real(Hx), real(Hy));
    grid on
    title(['Plot of front view for TM_', num2str(m), '_', num2str(n)]);
    legend('E-Field', 'H-Field');
    xlabel('x-dimension 0 to a');
    ylabel('y-dimension 0 to b=a/2');

    % Similar plots can be generated for Top View and Side View
end
