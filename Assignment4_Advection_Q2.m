clear all;
clc;

%% STEP-1: User-Input and Physical Properties Problem-2
L = 6.0; H = 1.0;           % Non-dimensional domain size 
imax = 42; jmax = 22;       % Grid points 
Beta = 1.2;                 % Grid stretching parameter
Re = 10; Pr = 1;            % Specified Re and Pr 
u = 1.0; v = 0.0;           % Slug flow velocity 
rho = 1.0; cp = 1.0;        
k = (rho * cp * u * H) / (Re * Pr); % Thermal Conductivity calculation
epsilon_st = 0.0001;        % Convergence criteria 
Dt = 0.05;                  % Time step size

disp('SELECT THE ADVECTION SCHEME\n');
disp('1. FOU'); 
disp('2. QUICK');
scheme = input('ENTER THE ADVECTION SCHEME (1 or 2): ');

%% STEP-2: Non-Uniform Grid Generation 
% X-direction: Finest near inlet, coarser towards outlet
xi_x = linspace(0, 1, imax-1);
num_x = (Beta+1) - (Beta-1)*((Beta+1)/(Beta-1)).^(1-xi_x);
den_x = ((Beta+1)/(Beta-1)).^(1-xi_x) + 1;
x = L * num_x ./ den_x;

% Y-direction: 
xi_y = linspace(0, 1, jmax-1);
sig_y = 2*xi_y - 1;
num_y = (Beta+1)*(((Beta+1)/(Beta-1)).^sig_y) - (Beta-1);
den_y = 2*(1 + ((Beta+1)/(Beta-1)).^sig_y);
y = H * num_y ./ den_y;

% Cell widths and centers calculation
Dx = diff(x); Dy = diff(y);
xc = [x(1), (x(1:end-1)+x(2:end))/2, x(end)];
yc = [y(1), (y(1:end-1)+y(2:end))/2, y(end)];
DX = [0, Dx, 0]; DY = [0, Dy, 0];

%% STEP-3: IC and BCs 
theta = zeros(imax, jmax);          % Initial condition θ=0 
theta(1, :) = 1.0;                  % Inlet boundary θ=1 
theta(:, 1) = 0.0; theta(:, jmax) = 0.0; % Wall boundaries θ=0 

%% STEP-4: Iterative Solver (Gauss-Seidel) 
unsteadiness = 1; n = 0;
while unsteadiness > epsilon_st
    n = n + 1;
    theta_old = theta;
    error_inner = 1;
    
    while error_inner > 0.0001
        theta_old_iter = theta;
        theta(imax, :) = theta(imax-1, :); % Outlet: Neumann BC 
        
        for i = 2:imax-1
            for j = 2:jmax-1
                % Diffusion (Conduction) terms
                De = k * DY(j) / (xc(i+1) - xc(i));
                Dw = k * DY(j) / (xc(i) - xc(i-1));
                Dn = k * DX(i) / (yc(j+1) - yc(j));
                Ds = k * DX(i) / (yc(j) - yc(j-1));
                
                % Advection terms (Upwind flux splitting)
                Ce = min(rho*u*DY(j), 0);
                Cw = max(rho*u*DY(j), 0);
                
                % Coefficient assembly
                ap0 = rho*DX(i)*DY(j)/Dt;
                ae = De - Ce;
                aw = Dw + Cw;
                an = Dn; as = Ds;
                ap = ap0 + ae + aw + an + as;
                
                % Deferred correction for QUICK Scheme 
                Qd = 0;
                if scheme == 2 && i > 2 && i < imax-1
                    curv = (theta_old_iter(i+1,j) - 2*theta_old_iter(i,j) + theta_old_iter(i-1,j));
                    Qd = 0.125 * rho * u * DY(j) * curv;
                end
                
                % Source term and GS Update 
                b = ap0 * theta_old(i,j) - Qd;
                theta(i,j) = (ae*theta(i+1,j) + aw*theta(i-1,j) + an*theta(i,j+1) + as*theta(i,j-1) + b) / ap;
            end
        end
        error_inner = max(max(abs(theta - theta_old_iter)));
    end
    unsteadiness = max(max(abs(theta - theta_old)))/Dt;
    fprintf('Time step: %d, Unsteadiness: %e\n', n, unsteadiness);
end

%% STEP-5: Plotting results 
% Figure 1: Contours 
figure(1); 
contourf(xc, yc, theta', 20); colorbar; colormap('jet');
title(['Steady State Temperature Contours (Scheme: ', num2str(scheme), ')']);
xlabel('X = x/H'); ylabel('Y');

% Figure 2: Profiles at axial locations 
figure(2);
hold on;
x_targets = [1, 2, 3, 5]; % Axial locations X=x/H 
markers = ['s', '*', 'o', 'd'];
for m = 1:4
    [~, idx] = min(abs(xc - x_targets(m)));
    plot(theta(idx, :), yc, ['-', markers(m)], 'DisplayName', ['X=', num2str(x_targets(m))]);
end
legend('Location', 'best'); grid on;
xlabel('\theta'); 
ylabel('Y'); 
title(['\theta(Y) Profiles at X=1,2,3,5 (Scheme: ', num2str(scheme), ')']);