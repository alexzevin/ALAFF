close all

% Set number of iterations to be performed
nk = 300

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
for i=1:N+2
    for j=1:N+2
        F( i,j ) = ...
            ( alpha^2 + beta^2 ) * pi^2 * sin( alpha * pi * x( i ) ) * sin( beta * pi * y( j ) );
    end
end

% Set the initial values at the mesh points
U = zeros( N+2, N+2 );

% set relationation parameter omega and initialize errors
omega_values = [1.0, 1.25, 1.5, 1.75, 1.9];
errors = zeros(length(omega_values), nk);

for w = 1:length(omega_values)  % loop through omega values
    omega = omega_values(w); % set omega to test
    U = zeros(N+2, N+2); % reset U

    for k = 1:nk
        U_old = U; % store old U
        
        % update all the interior points
        for i = 2:N+1
            for j = 2:N+1
                U_new = ( U( i, j-1 ) + U( i-1, j ) + U( i+1, j ) + U( i, j+1 ) + h^2 * F( i, j ) ) / 4;
                U(i,j) = (1 - omega) * U(i,j) + omega * U_new; % update with relaxation
            end
        end 
        errors(w, k) = max(max(abs(U - U_old))); % compute errors

        if mod(k, 10) == 0  
            mesh( x, y, U );
            title(['Iteration ', num2str(k), ', \omega = ', num2str(omega)]);
            axis( [ 0 1 0 1 -1.5 1.5 ]);
            drawnow;

            % wait to continue to the next iteration
            next = input( 'press RETURN to continue' );
        end
    end
end

% plot errors
figure; 
hold on; 
for w = 1:length(omega_values) 
    plot(1:nk, errors(w, :), 'DisplayName', ['\omega = ', num2str(omega_values(w))]);
end
hold off;
xlabel('Iteration');
ylabel('Error');
legend;
title('Error convergence for different \omega values');

