clc;
clear;
close all;

%% Lid-Driven Cavity Using Staggered Grid

% Thermophysical properties
rho = 1;
U   = 1;
cp  = 1;
V   = 0.5;    % only for stability estimate
L   = 1;
Re  = 100;
mu  = 1/Re;

imax = 32;
jmax = 32;

%% Step 2: Geometry
Dx = L / (imax - 2);
Dy = L / (jmax - 2);

dx = Dx;
dy = Dy;

Epsilon_st = 1e-6;

%% Staggered-grid data structures
u = zeros(imax-1, jmax);   % u-velocity at vertical faces
v = zeros(imax, jmax-1);   % v-velocity at horizontal faces
p = zeros(imax, jmax);     % pressure at cell centers

u_star    = zeros(imax-1, jmax);
v_star    = zeros(imax, jmax-1);
p_star_new = zeros(imax, jmax);

%% Step 3: Initial conditions / boundary conditions
u(:,1)    = 0;   % bottom wall
u(:,jmax) = U;   % top lid
v(1,:)    = 0;   % left wall
v(imax,:) = 0;   % right wall

u_star(:) = u(:);
v_star(:) = v(:);

Error_st = 10;

%% Step 4: Store previous step
u_prev_step = u;
v_prev_step = v;

%% Step 5: CFL Stability criteria
Dt_adv  = 1 / (abs(U)/Dx + abs(V)/Dy);
Dt_diff = 0.5 / (mu * (1/Dx^2 + 1/Dy^2));
Dt      = 0.6 * min(Dt_adv, Dt_diff);

%% Preallocate working arrays
% m_ux   = zeros(imax-1, jmax);
% m_ux_p = zeros(imax-1, jmax);
% m_ux_m = zeros(imax-1, jmax);
% a_ux   = zeros(imax-1, jmax);
% d_ux   = zeros(imax-1, jmax);
% 
% m_uy   = zeros(imax-1, jmax-1);
% m_uy_p = zeros(imax-1, jmax-1);
% m_uy_m = zeros(imax-1, jmax-1);
% a_uy   = zeros(imax-1, jmax-1);
% d_uy   = zeros(imax-1, jmax-1);
% 
% m_vx   = zeros(imax-1, jmax-1);
% m_vx_p = zeros(imax-1, jmax-1);
% m_vx_m = zeros(imax-1, jmax-1);
% a_vx   = zeros(imax-1, jmax-1);
% d_vx   = zeros(imax-1, jmax-1);
% 
% m_vy   = zeros(imax, jmax-1);
% m_vy_p = zeros(imax, jmax-1);
% m_vy_m = zeros(imax, jmax-1);
% a_vy   = zeros(imax, jmax-1);
% d_vy   = zeros(imax, jmax-1);
% 
% A_ux = zeros(imax-1, jmax);
% D_ux = zeros(imax-1, jmax);
% S_u  = zeros(imax-1, jmax);
% 
% A_vx = zeros(imax, jmax-1);
% D_vx = zeros(imax, jmax-1);
% S_v  = zeros(imax, jmax-1);

m_x_star = zeros(imax-1, jmax);
m_y_star = zeros(imax, jmax-1);

m_x_prime_new = zeros(imax-1, jmax);
m_y_prime_new = zeros(imax, jmax-1);

S_m_p_star = zeros(imax, jmax);
S_m_p_corr = zeros(imax, jmax);

%% Time marching
while Error_st >= Epsilon_st

    u = u_prev_step;
    v = v_prev_step;

    % Re-apply velocity boundary conditions
    u(:,1)    = 0;
    u(:,jmax) = U;
    v(1,:)    = 0;
    v(imax,:) = 0;

    % Initialize predicted fields from current step
    u_star = u;
    v_star = v;

    % Reset arrays every outer iteration
    m_x_star(:) = 0;
    m_y_star(:) = 0;
    S_m_p_star(:) = 0;
    S_m_p_corr(:) = 0;

    A_ux(:) = 0; D_ux(:) = 0; S_u(:) = 0;
    A_vx(:) = 0; D_vx(:) = 0; S_v(:) = 0;

    %% Step 6: Predictor step
    %% First: Predict the fluxes at all internal grid points

    % u-control volume, x-direction
    for i = 1:imax-2
        for j = 2:jmax-1
            m_ux(i,j)   = rho * (u(i+1,j) + u(i,j)) / 2;
            m_ux_p(i,j) = max(0, m_ux(i,j));
            m_ux_m(i,j) = min(0, m_ux(i,j));

            a_ux(i,j) = m_ux_p(i,j) * u(i,j) + m_ux_m(i,j) * u(i+1,j);
            d_ux(i,j) = mu * (u(i+1,j) - u(i,j)) / Dx;
        end
    end

    % u-control volume, y-direction
    for i = 2:imax-2
        for j = 1:jmax-1
            m_uy(i,j)   = rho * (v(i+1,j) + v(i,j)) / 2;
            m_uy_p(i,j) = max(0, m_uy(i,j));
            m_uy_m(i,j) = min(0, m_uy(i,j));

            a_uy(i,j) = m_uy_p(i,j) * u(i,j) + m_uy_m(i,j) * u(i,j+1);

            if j == 1 || j == jmax-1
                d_uy(i,j) = mu * (u(i,j+1) - u(i,j)) / (0.5 * Dy);
            else
                d_uy(i,j) = mu * (u(i,j+1) - u(i,j)) / Dy;
            end
        end
    end

    %% v-control volume, x-direction
    for i = 1:imax-1
        for j = 2:jmax-2
            m_vx(i,j)   = rho * (u(i,j+1) + u(i,j)) / 2;
            m_vx_p(i,j) = max(0, m_vx(i,j));
            m_vx_m(i,j) = min(0, m_vx(i,j));

            a_vx(i,j) = m_vx_p(i,j) * v(i,j) + m_vx_m(i,j) * v(i+1,j);

            if i == 1 || i == imax-1
                d_vx(i,j) = mu * (v(i+1,j) - v(i,j)) / (0.5 * Dx);
            else
                d_vx(i,j) = mu * (v(i+1,j) - v(i,j)) / Dx;
            end
        end
    end

    %% v-control volume, y-direction
    for i = 2:imax-1
        for j = 1:jmax-2
            m_vy(i,j)   = rho * (v(i,j+1) + v(i,j)) / 2;
            m_vy_p(i,j) = max(0, m_vy(i,j));
            m_vy_m(i,j) = min(0, m_vy(i,j));

            a_vy(i,j) = m_vy_p(i,j) * v(i,j) + m_vy_m(i,j) * v(i,j+1);
            d_vy(i,j) = mu * (v(i,j+1) - v(i,j)) / Dy;
        end
    end

    %% Second: Total advection, diffusion, and source terms

    % u-control volume
    for i = 2:imax-2
        for j = 2:jmax-1
            A_ux(i,j) = (a_ux(i,j) - a_ux(i-1,j)) * Dy + (a_uy(i,j) - a_uy(i,j-1)) * Dx;
            D_ux(i,j) = (d_ux(i,j) - d_ux(i-1,j)) * Dy + (d_uy(i,j) - d_uy(i,j-1)) * Dx;
            S_u(i,j)  = (p(i,j) - p(i+1,j)) * Dy;

            u_star(i,j)   = u(i,j) + (Dt / (rho * Dx * Dy)) * (D_ux(i,j) - A_ux(i,j) + S_u(i,j));
            m_x_star(i,j) = rho * u_star(i,j);
        end
    end

    % v-control volume
    for i = 2:imax-1
        for j = 2:jmax-2
            A_vx(i,j) = (a_vx(i,j) - a_vx(i-1,j)) * Dy + (a_vy(i,j) - a_vy(i,j-1)) * Dx;
            D_vx(i,j) = (d_vx(i,j) - d_vx(i-1,j)) * Dy + (d_vy(i,j) - d_vy(i,j-1)) * Dx;
            S_v(i,j)  = (p(i,j) - p(i,j+1)) * Dx;

            v_star(i,j)   = v(i,j) + (Dt / (rho * Dx * Dy)) * (D_vx(i,j) - A_vx(i,j) + S_v(i,j));
            m_y_star(i,j) = rho * v_star(i,j);
        end
    end

    % Boundary mass fluxes
    for i = 2:imax-1
        m_y_star(i,1)      = rho * v(i,1);
        m_y_star(i,jmax-1) = rho * v(i,jmax-1);
    end
    for j = 2:jmax-1
        m_x_star(1,j)      = rho * u(1,j);
        m_x_star(imax-1,j) = rho * u(imax-1,j);
    end

    %% Predicted mass source at internal grid points
    for i = 2:imax-1
        for j = 2:jmax-1
            S_m_p_star(i,j) = (m_x_star(i,j) - m_x_star(i-1,j)) * Dy ...
                            + (m_y_star(i,j) - m_y_star(i,j-1)) * Dx;
        end
    end

  %% Check predicted mass conservation
Smp_star = max(abs(S_m_p_star(:)));

if Smp_star < Epsilon_st
    u = u_star;
    v = v_star;
    p_star_new = p;
else
    %% Repeat correction until mass conservation is satisfied
    Epsilon_mass = 1e-6;
    Residual_Mass = Smp_star;
    maxCorrIter = 600;
    corrIter = 0;

    while Residual_Mass > Epsilon_mass && corrIter < maxCorrIter

        corrIter = corrIter + 1;

        %% Pressure / mass-source correction field
        p_prime_new = zeros(imax, jmax);
        p_prime_old = zeros(imax, jmax);

        Epsilon = 1e-6;
        error   = 10;

        aE = Dt * Dy / dx;
        aW = Dt * Dy / dx;
        aN = Dt * Dx / dy;
        aS = Dt * Dx / dy;
        aP = aE + aW + aN + aS;

        %% Solve correction equation using current mass imbalance
        while error > Epsilon
            p_prime_old = p_prime_new;

            for i = 2:imax-1
                for j = 2:jmax-1
                    p_prime_new(i,j) = ...
                        ( aE * p_prime_new(i+1,j) ...
                        + aW * p_prime_new(i-1,j) ...
                        + aN * p_prime_new(i,j+1) ...
                        + aS * p_prime_new(i,j-1) ...
                        - S_m_p_star(i,j) ) / aP;
                end
            end

            % Pressure-correction BCs
            p_prime_new(1,:)    = p_prime_new(2,:);
            p_prime_new(imax,:) = p_prime_new(imax-1,:);
            p_prime_new(:,1)    = p_prime_new(:,2);
            p_prime_new(:,jmax) = p_prime_new(:,jmax-1);

            % Reference point
            %p_prime_new(2,2) = 0;

            error = max(abs(p_prime_new(:) - p_prime_old(:)));
        end

        %% Flux correction
        m_x_prime_new(:) = 0;
        m_y_prime_new(:) = 0;

        for i = 2:imax-2
            for j = 2:jmax-1
                m_x_prime_new(i,j) = -Dt * (p_prime_new(i+1,j) - p_prime_new(i,j)) / dx;
            end
        end

        for i = 2:imax-1
            for j = 2:jmax-2
                m_y_prime_new(i,j) = -Dt * (p_prime_new(i,j+1) - p_prime_new(i,j)) / dy;
            end
        end

        %% Correct mass fluxes
        for i = 2:imax-2
            for j = 2:jmax-1
                m_x_star(i,j) = m_x_star(i,j) + m_x_prime_new(i,j);
            end
        end

        for i = 2:imax-1
            for j = 2:jmax-2
                m_y_star(i,j) = m_y_star(i,j) + m_y_prime_new(i,j);
            end
        end

        % Re-apply boundary mass fluxes
        for i = 2:imax-1
            m_y_star(i,1)      = rho * v(i,1);
            m_y_star(i,jmax-1) = rho * v(i,jmax-1);
        end
        for j = 2:jmax-1
            m_x_star(1,j)      = rho * u(1,j);
            m_x_star(imax-1,j) = rho * u(imax-1,j);
        end

        %% Update velocities from corrected mass fluxes
        for i = 2:imax-2
            for j = 2:jmax-1
                u_star(i,j) = m_x_star(i,j) / rho;
            end
        end

        for i = 2:imax-1
            for j = 2:jmax-2
                v_star(i,j) = m_y_star(i,j) / rho;
            end
        end

        % Re-apply velocity BCs
        u_star(:,1)    = 0;
        u_star(:,jmax) = U;
        v_star(1,:)    = 0;
        v_star(imax,:) = 0;

        %% Recompute mass imbalance after correction
        S_m_p_corr(:) = 0;

        for i = 2:imax-1
            for j = 2:jmax-1
                S_m_p_corr(i,j) = (m_x_star(i,j) - m_x_star(i-1,j)) * Dy ...
                                + (m_y_star(i,j) - m_y_star(i,j-1)) * Dx;
            end
        end

        Residual_Mass = max(abs(S_m_p_corr(:)));

        %% Update pressure incrementally
        p = p + p_prime_new;

        p(1,:)    = p(2,:);
        p(imax,:) = p(imax-1,:);
        p(:,1)    = p(:,2);
        p(:,jmax) = p(:,jmax-1);
        p(2,2)    = 0;

        %% Use corrected residual as next source for another correction if needed
        S_m_p_star = S_m_p_corr;

        fprintf('Correction iter = %d, Residual_Mass = %.6e\n', corrIter, Residual_Mass);
    end

    p_star_new = p;

    if Residual_Mass <= Epsilon_mass
        fprintf('Mass Conservation Satisfied\n');
    else
        fprintf('Warning: Mass conservation not fully satisfied after maxCorrIter\n');
    end
end

    %% Step 13: Update for next time step
    p = p_star_new;
    u = u_star;
    v = v_star;

    % Re-apply velocity boundary conditions
    u(:,1)    = 0;
    u(:,jmax) = U;
    v(1,:)    = 0;
    v(imax,:) = 0;

    % Pressure reference
    p(2,2) = 0;

    %% Step 14: Stopping criterion
    U_error = max(abs(u(:) - u_prev_step(:)));
    V_error = max(abs(v(:) - v_prev_step(:)));

    steady_error = max([U_error, V_error]);
    Error_st     = steady_error;

    if steady_error <= Epsilon_st
        fprintf('Solution Reached Steady State!\n');
        break;
    end

    u_prev_step = u;
    v_prev_step = v;
end

%% ============================================================
%% Post-processing
%% ============================================================

%% Cell-center velocity for plotting
xp = linspace(0, L, imax);
yp = linspace(0, L, jmax);
[Xp, Yp] = meshgrid(xp, yp);

Uc = zeros(imax, jmax);
Vc = zeros(imax, jmax);

% u at cell centers
for i = 2:imax-1
    for j = 1:jmax
        Uc(i,j) = 0.5 * (u(i,j) + u(i-1,j));
    end
end

% v at cell centers
for i = 1:imax
    for j = 2:jmax-1
        Vc(i,j) = 0.5 * (v(i,j) + v(i,j-1));
    end
end

% Boundary values at cell centers
Uc(1,:)    = 0;
Uc(imax,:) = 0;
Vc(:,1)    = 0;
Vc(:,jmax) = 0;

Umag  = sqrt(Uc.^2 + Vc.^2);
p_plot = p - p(2,2);

%% Velocity vector plot
figure
quiver(Xp', Yp', Uc, Vc, 'LineWidth', 1.0)
axis equal tight
xlabel('x')
ylabel('y')
title('Velocity Field')
figure
streamslice(Uc',Vc',1.1 )

axis equal tight
xlabel('x')
ylabel('y')
title('Streamline Plot')

%% Velocity magnitude contour
figure
contourf(Xp', Yp', Umag, 30, 'LineColor', 'none')
colorbar
colormap(jet)
axis equal tight
xlabel('x')
ylabel('y')
title('Velocity Magnitude Contour')

%% u-velocity contour
figure
contourf(Xp', Yp', Uc, 30, 'LineColor', 'none')
colorbar
colormap(jet)
axis equal tight
xlabel('x')
ylabel('y')
title('u-Velocity Contour')

%% v-velocity contour
figure
contourf(Xp', Yp', Vc, 30, 'LineColor', 'none')
colorbar
colormap(jet)
axis equal tight
xlabel('x')
ylabel('y')
title('v-Velocity Contour')

%% Pressure contour
figure
contourf(Xp', Yp', p_plot, 50, 'LineColor', 'none')
colorbar
colormap(jet)
axis equal tight
xlabel('x')
ylabel('y')
title('Pressure Contour')



%% Validation: u-velocity along x = 0.5
ghia_u_data = [
    1.0000  1.00000
    0.9766  0.84123
    0.9688  0.78871
    0.9609  0.73722
    0.9531  0.68717
    0.8516  0.23151
    0.7344  0.00332
    0.6172 -0.13641
    0.5000 -0.20581
    0.4531 -0.21090
    0.2813 -0.15662
    0.1719 -0.10150
    0.1016 -0.06434
    0.0703 -0.04775
    0.0625 -0.04192
    0.0547 -0.03717
    0.0000  0.00000
];

[~, mid_i] = min(abs(xp - 0.5));
u_center = Uc(mid_i, :);

figure
plot(u_center, yp, 'b-o', 'LineWidth', 1.2, 'DisplayName', 'My Solver')
hold on
plot(ghia_u_data(:,2), ghia_u_data(:,1), 'rs', ...
    'MarkerFaceColor', 'r', 'DisplayName', 'Ghia et al. (1982)')
grid on
xlabel('u-velocity')
ylabel('y')
title(['Validation at Re = ', num2str(Re), ': u at x = 0.5'])
legend('Location', 'best')

%% Validation: v-velocity along y = 0.5
ghia_v_data = [
    1.0000  0.00000
    0.9688 -0.05906
    0.9609 -0.07391
    0.9531 -0.08864
    0.9453 -0.10313
    0.9063 -0.16914
    0.8594 -0.22445
    0.8047 -0.24533
    0.5000  0.05454
    0.2344  0.17527
    0.2266  0.17507
    0.1563  0.16077
    0.0938  0.12317
    0.0781  0.10890
    0.0703  0.10091
    0.0625  0.09233
    0.0000  0.00000
];

[~, mid_j] = min(abs(yp - 0.5));
v_center = Vc(:, mid_j);

figure
plot(xp, v_center, 'b-o', 'LineWidth', 1.2, 'DisplayName', 'My Solver')
hold on
plot(ghia_v_data(:,1), ghia_v_data(:,2), 'rs', ...
    'MarkerFaceColor', 'r', 'DisplayName', 'Ghia et al. (1982)')
grid on
xlabel('x')
ylabel('v-velocity')
title(['Validation at Re = ', num2str(Re), ': v at y = 0.5'])
legend('Location', 'best')