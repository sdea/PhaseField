%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  --- SBM NO CONTACT ANGLE 3D     -------------------      %%%%%%%%%%%%% 
% REAL SIMULATIONS WITH SMOOTH BOUNDARY CONDTIONS
% Needed smoothed geometries

clc;
clear;
close all;


%% Set parameter simulation
h = 0.2; 
dt = 1e-6; 
T = 100;
begin = true;
smoothed = false;
SMobility = true;
BC = true;
test = false;
manual = true;
fromFunc = false;

% Paramater about mobility, gradient energy and energy funcitonal
W =  1;
Ms = 1;
eps = 1;
k = eps^2;

%% Form here the order parameter and SBM function from previous smoothing
Ctmp = load('Cini_3D.mat');
Ctmp = Ctmp.Cini_3D;

C = zeros(size(Ctmp));
C(Ctmp == 3) = 1;

Psi = load('phiYSZ.mat');
Psi = Psi.phiYSZ;


% Get dimensions
dim = size(C);
N = dim(1);

ivec = 3:N-2;

%% Compute dPsi/Psi to use afterwards in simulation
G = 1./(Psi+1);
[dPsix, dPsiy, dPsiz] = grad_3D(Psi,h);
inv_dPsix = G(ivec,ivec,ivec).*dPsix;
inv_dPsiy = G(ivec,ivec,ivec).*dPsiy;
inv_dPsiz = G(ivec,ivec,ivec).*dPsiz;

%% Filling to prevent mass losses
%  C(Psi > 0.5 & Psi < 0.95 & C > 0.05 & C < 0.95) = 1;




%% Case of no contact angle
M = ones(size(C));
dM = zeros(size(C));
c =1;
figure,
for i=1: round(T/dt)
    
    % Disp time step
    t = i*dt;
    disp(i);
    
    % Energy function and its derivative to respect C
    F = (W/2)*C.^2.*(1 - C).^2;
    dF = W*(1-C).*(1-2*C).*C;
    
    % Mobility function of C with boxcar function
    if (i > 500)
        if (SMobility)
            M = Ms*C.^2.*(1 - C).^2;
            dM = 2*Ms*C.*(1-C).*(1-2*C);
        end
    end
    
    % Intermediate terms
    [dCx,dCy,dCz] = grad_3D(C,h);
    [dFx,dFy,dFz] = grad_3D(dF,h);
    
    % Laplacian and derivative
    d2C = d2x_3D(C,h);
    
    
    % g = grad(nabla^2C)
    [dgx,dgy,dgz] = grad_3D(d2C,h);
    
    % Compute inner paranthesis
    Hx = dFx(3:end-2,3:end-2,3:end-2) - k*dgx;
    Hy = dFy(3:end-2,3:end-2,3:end-2) - k*dgy;
    Hz = dFz(3:end-2,3:end-2,3:end-2) - k*dgz;
    
    % Final computation
    t1 = M(ivec,ivec,ivec).*(d2C - d4x_3D(C,h));
    
    t2 = dM(5:end-4,5:end-4,5:end-4).*dCx(3:end-2,3:end-2,3:end-2).*Hx + ...
        dM(5:end-4,5:end-4,5:end-4).*dCy(3:end-2,3:end-2,3:end-2).*Hy + ...
        dM(5:end-4,5:end-4,5:end-4).*dCz(3:end-2,3:end-2,3:end-2).*Hz;
    
    t3 = M(5:end-4,5:end-4,5:end-4).*(inv_dPsix(3:end-2,3:end-2,3:end-2).*Hx + ...
                inv_dPsiy(3:end-2,3:end-2,3:end-2).*Hy + ... 
                inv_dPsiz(3:end-2,3:end-2,3:end-2).*Hz);
    
    
    
    % Time increment
    dCdt = t1(3:end-2,3:end-2,3:end-2) + t2 + t3;
    
    % Forward euler
    C(5:end-4,5:end-4,5:end-4) = C(5:end-4,5:end-4,5:end-4) + dt*dCdt;
    
 
    
    % Check on the internal energy functional
    if (mod(i,100) == 0)
        
        c = c+1;
        F_func = F(3:end-2,3:end-2,3:end-2) + (k/2)*dx_3D(C,h);
        F_func = sum(F_func(:));
        Ft(c) = F_func*h;
        imagesc(C(5:end-4,5:end-4,30));title(sprintf('F=%e, t=%e',F_func,t)); colormap(jet); drawnow();
    end
    
    if (mod(i,2000) == 0)
        filename = strcat('OutputMobility/C_',int2str(i));
        disp('Saving...');
        save(filename,'C');
    end
    
end
% 
% %% Algorithm for retrieving the phases
% tmpGeom = zeros(size(C));
% tmpGeom(C > 0.5) = 2;
% tmpGeom(Psi < 0.5) = 1;
% 
% %% Produce stl files
% padgeom = padarray(tmpGeom, [1 1 1], 0);
% 
% 
% 
% % Nickel
% Ni = zeros(size(padgeom));
% Ni(padgeom == 2) = 1;
% fv_ni = isosurface(Ni, 0.5);
% stlwrite('Ni_ini.stl', fv_ni);
% 
% % YSZ
% ysz = zeros(size(padgeom));
% ysz(padgeom == 1) = 1;
% fv_ysz = isosurface(ysz, 0.5);
% stlwrite('YSZ_ini.stl',fv_ysz);
% 

%% 
files = dir();

S = [files.datenum];
[~,S] = sort(S);
sfiles = {files(S).name};

%%
for i=3:2:36
    
    name = sfiles{i};
    disp(name)
    load(name, 'C');
    
    % Tmp geometry
    tmpGeom = zeros(size(C));
    tmpGeom(C > 0.5) = 2;
    tmpGeom(Psi < 0.5) = 1;
    
    % Pad geometry
    padgeom = padarray(tmpGeom(5:end-4,5:end-4,5:end-4), [1 1 1], 0);
    
    % % Nickel
    Ni = zeros(size(padgeom));
    Ni(padgeom == 2) = 1;
    fv_ni = isosurface(Ni, 0.5);
    stlwrite(strcat('Stl/Ni',num2str(i,'%04d'),'.stl'), fv_ni);
    
end
 
%% YSZ
ysz = zeros(size(padgeom));
ysz(padgeom == 1) = 1;
fv_ysz = isosurface(ysz, 0.5);
stlwrite('Stl/YSZ_ini.stl',fv_ysz);










