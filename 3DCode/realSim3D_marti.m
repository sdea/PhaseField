clc;
clear all;

addpath('M:\MATLAB\functions\Phase_field');

%% Init domain
load('test.mat');
geom = test;
slice_flipper_3D(geom);
hx = 0.2;   % Grid spacing in space (x direction)
hy = 0.2;   % Grid spacing in space (y direction)
hz = 0.2;   % Grid spacing in space (z direction)
dt = 1E-6;  % Time resolution
T = 1000;   % Total length time (n time points)
W = 1;      % Depth energy in Landau potential 
bcType = 1;  % Periodic boundary condtions
eps = 1;

% Trig surface mobility instead of bulk mobility
SMobility = false;

% Geometry flags
ni_pore = false;
psi_ysz = true;
smooth = false;
dmap = true;

% Contact angle
theta = (100/180)*pi;
cosTheta = cos(theta);

if (dmap)
    
    % Make here C function based on distance map
    % Distance map for YSZ-(ni and pore)
    yszNi = zeros(size(geom));
    yszNi(geom == 2) = 1;
    yszNi(geom == 1) = 1;
%     slice_flipper_3D(yszNi);

    yszPore = zeros(size(geom));
    yszPore(geom == 3) = 1;
    yszPore(geom == 2) = 1;   
%     slice_flipper_3D(yszPore);

    % Distance maps
    yszNi_dm = bwdist(~yszNi);
    yszPore_dm = bwdist(~yszPore);

    % Assign YSZ voxels to the respective phases based on their distance
    mask = zeros(size(geom));
    mask(geom == 2) = 1;

    tmp = geom;
    
    tmp(mask == 1 & yszPore_dm < yszNi_dm) = 1;
    tmp(mask == 1 & yszPore_dm >= yszNi_dm) = 3;
    
%     slice_flipper_3D(tmp);
    
    C = zeros(size(geom));
    C(tmp == 1) = 0;
    C(tmp == 3) = 1;
    
%     slice_flipper_3D(C);
end

YSZ = zeros(size(test));
YSZ(test == 2) = 1;

if (psi_ysz)
    dis_ysz = SDF(YSZ);
    delta = 4;
    Psi_ysz = double((1 - tanh((dis_ysz)/delta))/2);  
end

%% Get dimensions
dim = size(C);
N = dim(1);

ivec = 3:N-2;

%% Compute dPsi/Psi to use afterwards in simulation
Psi = Psi_ysz;
InvPsi = 1./(Psi+1);
[dxPsi, dyPsi, dzPsi] = grad3D(Psi,hx,hy,hz);
dxPsiOvPsi = InvPsi(ivec,ivec,ivec).*dxPsi;
dyPsiOvPsi = InvPsi(ivec,ivec,ivec).*dyPsi;
dzPsiOvPsi = InvPsi(ivec,ivec,ivec).*dzPsi;
modGradPsi = sqrt(dxPsi.^2 + dyPsi.^2 + dzPsi.^2);

%% Mobility function

% Constant mobility 
M = ones(size(C));
dM = zeros(size(C));
c =1;
figure,

for i = 1:round(T/dt);
    % Disp time step
    t = i*dt;
    disp(i);
    
    % Energy function and its derivative to respect C
    F = (W/2)*C.^2.*(1 - C).^2;
    dF = W*(1-C).*(1-2*C).*C;
    
    [dxM, dyM, dzM] = grad3D(M, hx, hy, hz);
    
    [dxC,dyC,dzC] = grad3D(C,hx,hy,hz);

    %Inner part of the equation, called g here 
    g = (dxPsiOvPsi.*dxC + dyPsiOvPsi.*dyC + dzPsiOvPsi.*dzC) +...
        d2xyz(C, hx, hy, hz) + ...
        (modGradPsi.*InvPsi(3:end-2,3:end-2,3:end-2).*...
        sqrt(2*F(3:end-2,3:end-2,3:end-2))./eps)*cosTheta;
    
    % Auxiliary term
    k = dF(3:end-2,3:end-2,3:end-2) - eps^2*g;
    
    % Finish here the equation
    [dxk, dyk, dzk] = grad3D(k, hx, hy, hz);
    dCdt =  M(5:end-4,5:end-4,5:end-4).*(dxPsiOvPsi(3:end-2,3:end-2,3:end-2).*dxk +...
            dyPsiOvPsi(3:end-2,3:end-2,3:end-2).*dyk + ...
            dzPsiOvPsi(3:end-2,3:end-2,3:end-2).*dzk) + ...
            M(5:end-4,5:end-4,5:end-4).*d2xyz(k, hx, hy, hz) + ...
            (dxM(3:end-2,3:end-2,3:end-2).*dxk + dyM(3:end-2,3:end-2,3:end-2).*dyk +...
            dzM(3:end-2,3:end-2,3:end-2).*dzk);
 
    % Forward euler integration
    C(5:end-4,5:end-4,5:end-4) = C(5:end-4,5:end-4,5:end-4) + dt*dCdt;
    
       if (bcType == 1)
        C(6:8,5:end-4,5:end-4) = C(2:4,5:end-4,5:end-4);
        C(end-3:end-1,5:end-4,5:end-4) = C(end-7:end-5,5:end-4,5:end-4);
        C(5:end-4,6:8,5:end-4) = C(5:end-4,2:4,5:end-4);
        C(5:end-4,end-3:end-1,5:end-4) = C(5:end-4,end-7:end-5,5:end-4);
        C(5:end-4,5:end-4,6:8) = C(5:end-4,5:end-4,2:4);
        C(5:end-4,5:end-4,end-3:end-1) = C(5:end-4,5:end-4,end-7:end-5);  
       end
   
    % Check on the internal energy functional
        if (mod(i,100) == 0)
            c = c+1;
            modGradC = sqrt(dxC.^2 + dyC.^2);
            F_func = F(3:end-2,3:end-2,3:end-2) + ((eps^2)/2)*modGradC.^2;
            F_func = sum(F_func(:));
            Ft(c) = F_func;
            imagesc(C(5:end-4,5:end-4,20));title(sprintf('F=%e, t=%e',F_func,t));colormap(jet);drawnow();
        end
   
       if (mod(i,1000) == 0)
               tmpGeom = zeros(size(C));
               tmpGeom(C > 0.5) = 2;
               tmpGeom(Psi < 0.5) = 1;
               
               % Pad geometry
               padgeom = padarray(tmpGeom(5:end-4,5:end-4,5:end-4), [1 1 1], 0);
    
               % % Nickel
               Ni = zeros(size(padgeom));
               Ni(padgeom == 2) = 1;
               fv_ni = isosurface(Ni, 0.5);
               stlwrite(strcat('Output02/Stl/Ni',num2str(i,'%04d'),'.stl'),fv_ni);
                
               filename = strcat('Output02/C_',int2str(i));
               disp('Saving...');
               save(filename,'C');
       end

end

%% YSZ
ysz = zeros(size(padgeom));
ysz(padgeom == 1) = 1;
fv_ysz = isosurface(ysz, 0.5);
stlwrite('Stl/YSZ_ini.stl',fv_ysz);