close all

% Set number of iterations to be performed
nk = 300;

% Set parameters alpha and beta
alpha = 2;
beta  = 3;

% Set the number of meshpoints so the interior has N x N such points
N = 50;

% Compute the distance between mesh points, in each direction
h = 1/(N+1);

% We will have arrays that capture the boundary as well as the interior
% meshpoints.  As a result, we need those arrays to be of size (N+2) x
% (N+2) so that those indexed 2:N+1, 2:N+1 represent the interior. 

% Compute the x-values at each point i,j, including the boundary
x = h * [ 0:N+1 ];   % Notice this creates a row vector

% Compute the y-values at each point i,j, including the boundary
y = h * [ 0:N+1 ];   % Notice this creates a row vector

% Create an array that captures the load at each point i,j
F = zeros(N+2, N+2);  % Initialize F array
for i=1:N+2
    for j=1:N+2
        F(i,j) = (alpha^2 + beta^2) * pi^2 * sin(alpha * pi * x(i)) * sin(beta * pi * y(j));
    end
end

% Set the initial values at the mesh points
U = zeros(N+2, N+2);

% Set the relaxation parameter omega (new line)
omega_values = [1.0, 1.25, 1.5, 1.75, 1.9];  % Different omega values for testing (new line)
errors = zeros(length(omega_values), nk);    % To store errors for each omega (new line)

% Storage for intermediate solution plots
solutions = cell(length(omega_values), 6);  % Store solutions at specific iterations

for w = 1:length(omega_values)  % Loop over different omega values (new line)
    omega = omega_values(w);    % Set current omega value (new line)
    U = zeros(N+2, N+2);  % Reset U for each omega (new line)

    for k = 1:nk
        U_old = U;  % Store old U for error computation (new line)
        
        % update all the interior points
        for i = 2:N+1
            for j = 2:N+1
                U_new = (U(i, j-1) + U(i-1, j) + U(i+1, j) + U(i, j+1) + h^2 * F(i, j)) / 4;  % Calculate new U value (new line)
                U(i,j) = (1 - omega) * U(i,j) + omega * U_new;  % Update U with relaxation
            end
        end 
        % Compute error
        errors(w, k) = max(max(abs(U - U_old)));  % Calculate max error (new line)
        
        % Store solutions at specific iterations for later plotting
        if ismember(k, [50, 100, 150, 200, 250, 300])
            solutions{w, k/50} = U;  % Store solution (new line)
        end
        
        % Plotting (optional)
        if mod(k, 50) == 0  % Plot every 50 iterations (new line)
            figure;
            mesh(x, y, U);
            title(['Iteration ', num2str(k), ', \omega = ', num2str(omega)]);  % Update title with current omega (new line)
            axis([0 1 0 1 -1.5 1.5]);
            saveas(gcf, ['Iteration_', num2str(k), '_omega_', num2str(omega), '.png']);  % Save figure (new line)
        end
    end
end

% Plot errors for different omega values (new section)
figure;  % Create new figure for error plots (new line)
hold on;  % Hold on to add multiple plots (new line)
for w = 1:length(omega_values)  % Loop over omega values (new line)
    plot(1:nk, errors(w, :), 'DisplayName', ['\omega = ', num2str(omega_values(w))]);  % Plot error for current omega (new line)
end
hold off;  % Release hold (new line)
xlabel('Iteration');  % Label x-axis (new line)
ylabel('Error');  % Label y-axis (new line)
legend;  % Add legend to plot (new line)
title('Error convergence for different \omega values');  % Add title to plot (new line)
saveas(gcf, 'Error_Convergence.png');  % Save error plot (new line)

% Plot intermediate solutions for each omega
for w = 1:length(omega_values)
    figure;
    for i = 1:6
        subplot(2, 3, i);
        mesh(x, y, solutions{w, i*50});
        title(['Iteration ', num2str(i*50), ', \omega = ', num2str(omega_values(w))]);
        axis([0 1 0 1 -1.5 1.5]);
    end
    saveas(gcf, ['Solutions_omega_', num2str(omega_values(w)), '.png']);  % Save solution plot (new line)
end
