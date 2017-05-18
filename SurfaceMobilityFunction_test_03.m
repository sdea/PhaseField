clc;
clear;

%% %% Set parameters

N = 50;    % Number of grid points 
hx = 0.2;   % Grid spacing in space x
hy = 0.2;   % Grid spacing in space y
dt = 1e-5;  % Time resolution
T = 1000;   % Total length time 
 
Q = 1/8;      % Depth energy in Landau potential 
bcType = 1; % Trig periodic boundary condition for spinoidal decomposition
eps = 1;

%% Geometry
geom = zeros(N+4,N+4);
geom((N+4)/2:end,:) = 1;
geom(1:(N+4)/2-1,(1:(N+4)/2)) = 2;
figure 
imagesc(geom) 
colormap(gray)

% Contact angle
theta = (90/180)*pi;
cosTheta = cos(theta);


%% Construct C function using the interface width
x = 1:N+4;
int_width = eps*sqrt(2/Q);
x0 = (N+4)/2;
C = repmat(((1-tanh((x-x0)/int_width))/2),N+4,1);
% figure, imagesc(C), colormap(jet); title C;


%% Construct Psi function'x = 1:N+4;
delta = 4;  
x0 = (N+4)/2;    
Psi = repmat(((1-tanh((x-x0)/delta))/2)',1,N+4);
% figure, imagesc(Psi), colormap(jet); title Psi;


%% Psi
% compute here the gradient of Psi and the term dPsi/Psi
invPsi = 1./Psi;
[dxPsi, dyPsi] = grad2D(Psi, hx, hy);
modGradPsi = sqrt(dxPsi.^2 + dyPsi.^2); 
dxPsiovPsi = invPsi(3:end-2,3:end-2).*dxPsi;
dyPsiovPsi = invPsi(3:end-2,3:end-2).*dyPsi;

%% Mobility function
dMC = zeros(N+4,N+4);

%% Solver

c = 1;
ca = 1;
MNi_P = 6;

for i=1:round(T/dt)    

    t = i*dt;   % Value of each istant
    disp(i);
       
      
    % Energy function and its derivative to respect C
    F = (Q/4)*C.^2.*(1 - C).^2;
    dF = Q/2*(1-C).*(1-2*C).*C;
    
    
    M = M_phi_psi(MNi_P, C, Psi);
    dMC = 2*C-4*C.^3;
%     dMC = 1 - 2*C;
    
    % Cut off (square box funciton)
    dMC(dMC > 0.9) = 0;
    dMC(dMC < 0.1) = 0;
    
    % 4th order mobility function
    % dMC = 2*C-4*C.^3;
    
    % if G is present (interpolation function)
    G = (Psi.^6).*(10*Psi.^2-15*Psi+6);
    dGPsi = (Psi.^5).*(80*Psi.^2-105*Psi+36);
    
%     G = ones(N+4,N+4);
%     dGPsi = zeros(N+4,N+4);
    
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
        (MNi_P*((dMC(5:end-4,5:end-4).*dxC(3:end-2,3:end-2).*G(5:end-4,5:end-4))+...
        (M(5:end-4,5:end-4).*dGPsi(5:end-4,5:end-4).*dxPsi(3:end-2,3:end-2))).*dxk + ...
        MNi_P*((dMC(5:end-4,5:end-4).*dyC(3:end-2,3:end-2).*G(5:end-4,5:end-4))+...
        (M(5:end-4,5:end-4).*dGPsi(5:end-4,5:end-4).*dyPsi(3:end-2,3:end-2))).*dyk);
          

        
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

if (mod(i,1000) == 0)
    ca = ca+1;
    % contact angle based on Chen et al.
    
    tpb_region = zeros(size(C(5:end-4,5:end-4)));
    tpb_region( C(5:end-4,5:end-4) > 0.1 &  C(5:end-4,5:end-4) < 0.9 ...
        & Psi(5:end-4,5:end-4) > 0.1 & Psi(5:end-4,5:end-4) < 0.9 ) = 1;
    den = (1./modGradC).*(1./modGradPsi);
    num = (dxC.*dxPsi + dyC.*dyPsi);
    M_contact_cos = num.*den;
    vec_contact = M_contact_cos(tpb_region == 1);
    contact_cos = mean(vec_contact);
    measured_theta_chen(ca) = acosd(-contact_cos);
end  
end

%% Perform some check on the contact angle
%  Indentifying region of TPB

tpb_region = zeros(size(C(5:end-4,5:end-4)));
tpb_region( C(5:end-4,5:end-4) > 0.3 &  C(5:end-4,5:end-4) < 0.7 ...
    & Psi(5:end-4,5:end-4) > 0.3 & Psi(5:end-4,5:end-4) < 0.7 ) = 1;

figure, imagesc(tpb_region), colormap(gray)

%% Contact angle based on Chen et al.

modGradC = sqrt(dxC.^2 + dyC.^2);
den = (1./modGradC).*(1./modGradPsi);
num = (dxC.*dxPsi + dyC.*dyPsi);
M = num.*den;
vec_contact = M(tpb_region == 1);
contact_cos = mean(vec_contact);

measured_theta = acosd(-contact_cos);
disp(measured_theta);

%% Retrieve phases
newgeom = zeros(size(C));
newgeom(C > 0.5) = 2;
newgeom(Psi < 0.5) = 1;
figure, imagesc(newgeom(5:end-4,5:end-4)), colormap(gray)




