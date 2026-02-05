%% FVM For Multi-Solid

Lx = 1;
Ly = 1;
imax = 20;
jmax = 20;
x=linspace(0,Lx,imax); 
y=linspace(0,Ly,jmax);
T = zeros(imax,jmax);
T0=30;
% Initialize temperature for the first time step
T(:,:) = T0;

T_left = 100;
T_top = 400;
T_right = 300;
T_bottom = 200;
Q_gen = 0;

Dtc = T_top-T_left ;
%Dt = 10;


T(1,:)=T_left;    
T(end,:)= T_right;   
T(:,end)=T_top;     
T(:,1)= T_bottom;  



Dt    = 5;
tEnd  = 10000;                 % total physical time
nStep = round(tEnd / Dt);

%% Meshing Infor

%% Centroid

xc(2:imax-1)=(x(2:imax-1)+x(1:imax-2))/2; %% Interior
xc(1)=x(1);
xc(imax)=x(imax-1);

yc(2:jmax-1) = (y(2:jmax-1)+y(1:jmax-2))/2;
yc(1)=y(1);
yc(jmax)=y(jmax-1);

%% Cell width and Height
for i=2:imax-1
    Dx(i)=x(i)-x(i-1);
end
for j=2:jmax-1
    Dy(j)=y(j)-y(j-1);
end

% Cell center Distance
  for i=1:imax-1
      dx(i)= (xc(i+1)-xc(i));
  end
  for j=1:jmax-1
      dy(j)=yc(j+1)-yc(j);
  end









%% Material ID Matrix
matID = zeros(imax,jmax);
for j = 1:jmax
    for i = 1:imax
        if x(i) <= 0.5*Lx && y(j) <= 0.5*Ly
            matID(i,j) = 1;
        elseif x(i) > 0.5*Lx && y(j) <= 0.5*Ly
            matID(i,j) = 2;
        elseif x(i) <= 0.5*Lx && y(j) > 0.5*Ly
            matID(i,j) = 3;
        else
            matID(i,j) = 4;
        end
    end
end
% Material properties 
rho_s = [8900, 7900, 8000, 2700];   
cp_s  = [385,  452,  502,  896];    
k_s   = [388,  72, 16.2,  220];    

%% Construct  Property Matrix
rho = zeros(imax,jmax);
cp  = zeros(imax,jmax);
k   = zeros(imax,jmax);
%% Assignment Property to Material
for j = 1:jmax
    for i = 1:imax
        id = matID(i,j);     % material number (1–4)
        rho(i,j) = rho_s(id);
        cp(i,j)  = cp_s(id);
        k(i,j)   = k_s(id);
    end
end

aE  = zeros(imax,jmax);
aW  = zeros(imax,jmax);
aN  = zeros(imax,jmax);
aS  = zeros(imax,jmax);
aP  = zeros(imax,jmax);
aP0 = zeros(imax,jmax);
b   = zeros(imax,jmax);
%% Compute LAE Coefficient
for j = 2:jmax-1
    for i = 2:imax-1

        aP0(i,j) = rho(i,j)*cp(i,j)*Dx(i)*Dy(j)/Dt;

        kE = 2*k(i,j)*k(i+1,j)/(k(i,j)+k(i+1,j));
        kW = 2*k(i,j)*k(i-1,j)/(k(i,j)+k(i-1,j));
        kN = 2*k(i,j)*k(i,j+1)/(k(i,j)+k(i,j+1));
        kS = 2*k(i,j)*k(i,j-1)/(k(i,j)+k(i,j-1));

        aE(i,j) = kE * Dy(j) / dx(i);
        aW(i,j) = kW * Dy(j) / dx(i-1);
        aN(i,j) = kN * Dx(i) / dy(j);
        aS(i,j) = kS * Dx(i) / dy(j-1);

        aP(i,j) = aP0(i,j) + aE(i,j) + aW(i,j) ...
                               + aN(i,j) + aS(i,j);
    end
end

 tol = 1e-4;
 iter    = 0;
 iterMax = 10000;
 unsteadiness_nd=1;
 epsilon_st=1e-4;
 alpha = zeros(imax,jmax);


for j = 2:jmax-1
    for i = 2:imax-1
        alpha(i,j) = k(i,j) / (rho(i,j)*cp(i,j));
    end
end
%% Now-Time Marching




%% Tecplot Output 

fid = fopen('D:/M.TECH/SEM2/CFD/Transient_MultiSolid.dat','w');

fprintf(fid,'TITLE = "Transient Multi-Solid Heat Conduction"\n');
fprintf(fid,'VARIABLES = "X" "Y" "T" "MatID"\n');

writeEvery = 5;

for n = 1:nStep
    time = n * Dt;
    iter = iter + 1;
   T(1,:)=T_left;    
  T(end,:)= T_right;   
  T(:,end)=T_top;     
  T(:,1)= T_bottom; 
    
    % Save previous time step
    T_old = T;
    % Right Hand Side
    for j = 2:jmax-1
        for i = 2:imax-1
            b(i,j) = aP0(i,j)*T_old(i,j) + Q_gen*Dx(i)*Dy(j);
        end
    end

%% Iteration for each time step using GS
Error = 1;
while Error>tol
   
  T(1,:)=T_left;    
  T(end,:)= T_right;   
  T(:,end)=T_top;     
  T(:,1)= T_bottom;
    T_old_iter = T; 
    for j = 2:jmax-1
        for i = 2:imax-1
            T(i,j) = ( ...
             aE(i,j) * T(i+1,j) + ...
            aW(i,j) * T(i-1,j) + ...
            aN(i,j) * T(i,j+1) + ...
            aS(i,j) * T(i,j-1) + ...
            b(i,j) ) / aP(i,j);
        end
    end
    Error = max(max(abs(T-T_old_iter)));
end
%% Write to Tecplot (TRANSIENT-CORRECT)
if mod(n, writeEvery) == 0

    fprintf(fid, ...
        'ZONE STRANDID=1, SOLUTIONTIME=%g, I=%d, J=%d, DATAPACKING=POINT\n', ...
        time, imax, jmax);

    for j = 1:jmax
        for i = 1:imax
            fprintf(fid,'%e %e %e %d\n', ...
                    x(i), y(j), T(i,j), matID(i,j));
        end
    end
end




%% Compute Steady State Criteria( This is Professor Atul Sharma Method)

%unsteadiness = max(max(abs(T - T_old))) / Dt;
%%alpha_min = min(alpha(alpha > 0));

%unsteadiness_nd = ...
 %%   * (Lx^2 / (alpha_min * Dt));
 
 




 %% This is generic method


 % --- Unsteadiness ---
    unsteadiness_nd = ...
        ( max(max(abs(T - T_old))) / Dtc ) ...
        * ( Lx^2 / (alpha_min * Dt) );

    % --- Monitor ---
    if mod(n,10) == 0
        fprintf('Step %5d | Time = %.2f | Unsteadiness = %.3e\n', ...
                n, time, unsteadiness_nd);
    end

    % --- Optional early stop ---
    if unsteadiness_nd < epsilon_st
        fprintf('Steady state reached at t = %.2f s\n', time);
        break
    end
 
 end


fclose(fid);

figure
contourf(x, y, T', 100, 'LineColor','none')
colorbar
colormap("jet")
xlabel('x')
ylabel('y')
title('Temperature distribution')
axis equal tight