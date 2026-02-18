%% STEP-1: User Input
clc; clear;

rho = 1000.0; 
cp = 4180.0;
Lx = 1.0;
Ly = 1.0;          
imax = 52;
jmax = 52;        % Grid Size
T0 =50;
T_wb =100; 
T_sb = 0.0; 

u = 1; 
v = 1; 
epsilon_st = 1e-5;  
Delta_Tc=100;
u_c=1;
disp('SELECT THE ADVECTION SCHEME (1/2/3)');
disp('1. FOU');
disp('2. SOU'); 
disp('3. QUICK');

scheme = input('ENTER THE ADVECTION SCHEME : ');

%% STEP-2: Geometrical Parameter & Time Step
Dx = Lx/(imax-2);
Dy = Ly/(jmax-2);
Dt = 1 / (abs(u)/Dx + abs(v)/Dy); % FOU scheme

if scheme == 2
    Dt = (2/4)*Dt; % SOU scheme
elseif scheme == 3
    Dt = (4/9)*Dt; % QUICK scheme
end

%% STEP-3: IC and BCs
T = zeros(imax, jmax);
T(2:imax-1, 2:jmax-1) = T0;
T(1, :) = T_wb; T(:, 1) = T_sb;

%% STEP-4: Calculate mass flux
mx = rho*u; my = rho*v;
mx_p = max(mx,0); mx_m = min(mx,0);
my_p = max(my,0); my_m = min(my,0);

%% FUNCTION: Compute face-center temperature
function [Tf_p, Tf_m] = Temp_f(w,T1,T2,T3,T4)
    Tf_p = w(1)*T3 + w(2)*T2 + w(3)*T1;
    Tf_m = w(1)*T2 + w(2)*T3 + w(3)*T4;
end

%% STEP-5: Time-Marching for Explicit Unsteady State
unsteadiness_nd = 1; n = 0;

while unsteadiness_nd >= epsilon_st
    n = n+1;
    
    % Non-Dirichlet BC: FOU at outlet
    T(imax, :) = T(imax-1, :);
    T(:, jmax) = T(:, jmax-1);
    
    T_old = T; % previous step
    
    %% Compute enthalpy flux in x-direction
    hx_old = zeros(imax-1,jmax-1);
    for j = 2:jmax-1
        for i = 1:imax-1
            if i==1 || i==imax-1
                Te_p = T_old(i,j); Te_m = T_old(i+1,j);
            else
                w = weights(scheme);
                [Te_p, Te_m] = Temp_f(w, T_old(i-1,j), T_old(i,j), T_old(i+1,j), T_old(i+2,j));
            end
            hx_old(i,j) = cp*(mx_p*Te_p + mx_m*Te_m);
        end
    end
    
    %% Compute enthalpy flux in y-direction
    hy_old = zeros(imax-1,jmax-1);
    for j = 1:jmax-1
        for i = 2:imax-1
            if j==1 || j==jmax-1
                Tn_p = T_old(i,j); Tn_m = T_old(i,j+1);
            else
                w = weights(scheme);
                [Tn_p, Tn_m] = Temp_f(w, T_old(i,j-1), T_old(i,j), T_old(i,j+1), T_old(i,j+2));
            end
            hy_old(i,j) = cp*(my_p*Tn_p + my_m*Tn_m);
        end
    end
    
    %% Update temperature
    for j = 2:jmax-1
        for i = 2:imax-1
            Q_adv_old = ((hx_old(i,j)-hx_old(i-1,j))*Dy) + ((hy_old(i,j)-hy_old(i,j-1))*Dx);
            T(i,j) = T_old(i,j) - (Dt/(rho*cp*Dx*Dy)) * Q_adv_old;
        end
    end
    
    %% Convergence criterion
    unsteadiness_nd = (Lx/(u_c*Delta_Tc))*max(max(abs(T - T_old))) / Dt;
    fprintf('Time step no. %5d, unsteadiness_nd = %8.4e\n', n, unsteadiness_nd);
end



%% STEP-6: Plotting the Temperature Field
x = linspace(0,Lx,imax);
y = linspace(0,Ly,jmax);
[X,Y] = meshgrid(x,y);

figure;
contourf(X,Y,T',60,'LineColor','none');   % transpose T to match meshgrid
colorbar;
colormap('jet');
xlabel('X'); 
ylabel('Y'); 
title('Temperature Contours');


%% Exact Solution and Centerline Comparison

%x = linspace(0,Lx,imax);
%y = linspace(0,Ly,jmax);

% Find index closest to x = 0.5
[~, i_center] = min(abs(x - 0.5));

T_center = T(i_center,:);

% Exact solution at x = 0.5 (REVERSED)
T_exact = zeros(size(y));
T_exact(y >= 0.5) = 100;
T_exact(y < 0.5)  = 0;

figure;
plot(y, T_center,'LineWidth',2,Marker='*');
hold on;
plot(y, T_exact,'k--','LineWidth',2);

xlabel('x');
ylabel('Temperature');
legend('Numerical','Exact');
title(['Vertical Centerline Profile (x=0.5), Scheme = ', num2str(scheme)]);
grid on;


%% Function Weight for Advection
function w = weights(scheme)
    switch scheme
        case 1
            w = [0,1,0];            % FOU
        case 2
            w = [0,1.5,-0.5];       % SOU
        case 3
            w = [3/8,6/8,-1/8];     % QUICK
        otherwise
            error('Invalid scheme selected.');
    end
end
