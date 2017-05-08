clc;
clear all;
addpath('M:\MATLAB\functions\Phase_field');

%% %% Set parameters

N = 100;    % Number of grid points 
hx = 0.2;   % Grid spacing in space x
hy = 0.2;   % Grid spacing in space y
dt = 1e-5;  % Time resolution
T = 1000;    % Total length time 
 
Q = 1;      % Depth energy in Landau potential 
bcType = 1; % Trig periodic boundary condition for spinoidal decomposition
eps = 1;

%% Geometry

geom = zeros(N+4,N+4);
geom((N+4)/2:end,:) = 1;
geom(1:(N+4)/2-1,(1:(N+4)/2)) = 2;
figure 
imagesc(geom) 
colormap(gray)

% If true use surface mobility instead of bulck mobility
SMobility = false;

% Geometry flags
ni_pore = false;
ni_ysz = true;
smooth = false;
dmap = true;

% Contact angle
theta = (100/180)*pi;
cosTheta = cos(theta);

%% Order parameter
if (ni_ysz)
    x = 1:N+4;
    delta = 4;%10;  %(sqrt(2/Q))*k;
    x0 = (N+4)/2;    
    Psi = repmat(((1-tanh((x-x0)/delta))/2)',1,N+4);
    
        if (smooth)
            C = zeros(N+4,N+4);
            x = 1:N+4;
            delta = 20;
            x0 = (N+4)/2;
                for i=1:length(x)
                    C(i,:) = (1 - tanh((x - x0)/delta))/2; 
                end
        end
    
end

if (dmap)
    
    % Make here C function based on distance map
    % Distance map for YSZ-(ni and pore)
    yszNi = zeros(size(geom));
    yszNi(geom == 1) = 1;
    yszNi(geom == 0) = 1;

    yszPore = zeros(size(geom));
    yszPore(geom == 1) = 1;
    yszPore(geom == 2) = 1;   

    % Distance maps
    yszNi_dm = bwdist(~yszNi);
    yszPore_dm = bwdist(~yszPore);

    % Assign YSZ voxels to the respective phases based on their distance
    mask = zeros(size(geom));
    mask(geom == 1) = 1;

    tmp = geom;
    
    tmp(mask == 1 & yszPore_dm < yszNi_dm) = 0;
    tmp(mask == 1 & yszPore_dm >= yszNi_dm) = 2;

    figure
    imagesc(tmp)
    colormap(gray);
    
    C = zeros(N+4,N+4);
    C(tmp == 0) = 0;
    C(tmp == 2) = 1;
end

%% Psi
% compute here the gradient of Psi and the term dPsi/Psi
invPsi = 1./Psi;
[dxPsi, dyPsi] = grad2D(Psi, hx, hy);
modGradPsi = sqrt(dxPsi.^2 + dyPsi.^2); 
dxPsiovPsi = invPsi(3:end-2,3:end-2).*dxPsi;
dyPsiovPsi = invPsi(3:end-2,3:end-2).*dyPsi;

%% Mobility function

% Constant mobility 
M = ones(N+4,N+4);
% M(M == 1) = 2;

% Surface mobility funtion
if (SMobility)
        M = M_phi_psi(1,C,Psi);
end
%% Solver

c = 1;
ca = 1;
old_meas_theta = 0;

for i=1:50000%round(T/dt)
    
    t = i*dt;   % Value of each istant
    %disp(i);
    
    % Energy function and its derivative to respect C
    F = (Q/4)*C.^2.*(1 - C).^2;
    dF = Q/2*(1-C).*(1-2*C).*C;
    
    [dxM, dyM] = grad2D(M, hx, hy);
    
    %Inner part of the equation, called g here 
    [dxC, dyC] = grad2D(C, hx, hy);
    g = (dxPsiovPsi.*dxC + dyPsiovPsi.*dyC) + d2xy(C, hx, hy) + ...
        (modGradPsi.*invPsi(3:end-2,3:end-2).*sqrt(2*F(3:end-2,3:end-2))./eps)*cosTheta;
    
    % Auxiliary term
    k = dF(3:end-2,3:end-2) - eps^2*g;
    
    % Finish here the equation
    [dxk, dyk] = grad2D(k, hx, hy);
    dCdt =  M(5:end-4,5:end-4).*(dxPsiovPsi(3:end-2,3:end-2).*dxk + dyPsiovPsi(3:end-2,3:end-2).*dyk) + ...
            M(5:end-4,5:end-4).*d2xy(k, hx, hy) + ...
            (dxM(3:end-2,3:end-2).*dxk + dyM(3:end-2,3:end-2).*dyk);
        
    % Forward euler integration
    C(5:end-4,5:end-4) = C(5:end-4,5:end-4) + dt*dCdt;
    
   if (bcType == 1)
        C(2:4,5:end-4) = C(6:8,5:end-4);
        C(end-3:end-1,5:end-4) = C(end-7:end-5,5:end-4);
        C(5:end-4, 6:8) = C(5:end-4,2:4);
        C(5:end-4, end-3:end-1) = C(5:end-4, end-7:end-5);
        
   end
    
       % Check on the internal energy functional
if (mod(i,1000) == 0)
        c = c+1;
        modGradC = sqrt(dxC.^2 + dyC.^2);
        F_func = F(3:end-2,3:end-2) + ((eps^2)/2)*modGradC.^2;
        F_func = sum(F_func(:));
        Ft(c) = F_func;
        imagesc(C(5:end-4,5:end-4));title(sprintf('F=%e, t=%e',F_func,t));colormap(jet); drawnow();
end

if (mod(i,100) == 0)
        ca = ca+1;
        % contact angle based on Chen et al.
        
        tpb_region = zeros(size(C(5:end-4,5:end-4)));
        tpb_region( C(5:end-4,5:end-4) > 0.1 &  C(5:end-4,5:end-4) < 0.9 ...
        & Psi(5:end-4,5:end-4) > 0.1 & Psi(5:end-4,5:end-4) < 0.9 ) = 1;
        modGradC = sqrt(dxC.^2 + dyC.^2);
        den = (1./modGradC).*(1./modGradPsi);
        num = (dxC.*dxPsi + dyC.*dyPsi);
        M_contact_cos = num.*den;
        vec_contact = M_contact_cos(tpb_region == 1);
        contact_cos = mean(vec_contact);
        measured_theta_chen(ca) = acosd(-contact_cos);
        
        new_meas_theta = measured_theta_chen(ca);
        disp (new_meas_theta)
        
        error(ca) = abs((new_meas_theta-old_meas_theta)/new_meas_theta);
        %disp(error);
        old_meas_theta = new_meas_theta;
        
    end
end

% %% Contact angle based on Chen et al.
% 
% modGradC = sqrt(dxC.^2 + dyC.^2);
% den = (1./modGradC).*(1./modGradPsi);
% num = (dxC.*dxPsi + dyC.*dyPsi);
% M = num.*den;
% vec_contact = M(tpb_region == 1);
% contact_cos = mean(vec_contact);
% 
% measured_theta = acosd(-contact_cos);
% disp(measured_theta);
% 
% %% Retrieve phases
% newgeom = zeros(size(C));
% newgeom(C > 0.5) = 2;
% newgeom(Psi < 0.5) = 1;
% figure, imagesc(newgeom(5:end-4,5:end-4)), colormap(gray)


