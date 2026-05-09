%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  AE 140
%      Aerospace Rigid Body Dynamics
% SpaceX Starship Earth to Earth Simultaion
%   Samuel Thomas Joseph, Joseph Joseph
%              May 11, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Constants and Variables
m_r = 5270000;                   % Total Mass of Rocket (kg)
g = 9.81;                        % Gravitational Constant (m/s^2)
r_earth = 6371000;               % Radius of Earth (m)
L = 55.74;                       % Distance from Base of Rocket to Center of Mass (m)
c_p = 55.09;                     % Location of Center of Pressure (m)
CD = 0.6;                        % Coefficient of Drag
A = 63.61;                       % Frontal Nose Area for Drag (m^2)
D_const = (0.5) * CD * rho * A;  % Drag Constant 

%% Inertia Dyadics (Tensors)
% Cone Dimension and Mass Properties 
m_cone = 397880;       % Mass of Cone (kg)
h_cone = 13.71;        % Height of Cone(m)
r_cone = 4.50;         % Radius of Cone (m)
z_cone = 50.87;        % Transfer Distance from Reference Point to COM of Cone (m)

% Cylinder Dimension and Mass Properties 
m_cylinder = 1052120;  % Mass of Cylinder (kg)
h_cylidner = 103.18;   % Height of Cylinder (m)
r_cylinder = 4.50;     % Radius of Cylinder (m)
z_cylinder = 4.15;     % Transfer Distance from Reference Point to COM of Cylinder (m)

% Cone Moments of Inertia
Ixx_Co = m_cone * ((3/20) * h_cone^2 + (3/80) * r_cone^2) + m_cone * z_cone^2;
Iyy_Co = m_cone * ((3/20) * h_cone^2 + (3/80) * r_cone^2) + m_cone * z_cone^2;
Izz_Co = m_cone * (3/10) * r_cone^2;

% Cylinder Moments of Inertia
Ixx_Cy = (1/12) * m_cylinder * (h_cylidner^2 + 3 * r_cylinder^2) + m_cylinder * z_cylinder^2;
Iyy_Cy = (1/12) * m_cylinder * (h_cylidner^2 + 3 * r_cylinder^2) + m_cylinder * z_cylinder^2;
Izz_Cy = (1/2) * m_cylinder * r_cylinder^2;

% Moments of Inertia Matrix
Icone = [Ixx_Co      0         0;
         0         Iyy_Co      0;
         0           0       Izz_Co];

Icylinder = [Ixx_Cy      0         0;
             0         Iyy_Cy      0;
             0           0       Izz_Cy];

% Inertia Dyadic of Rocket 
I = Icone + Icylinder;
Ixx = I(1,1);
Iyy = I(2,2);
Izz = I(3,3);

%% Atmospheric Model for Density 

% Altitude in terms of Cartesian Coordinates
h = (x^2 + y^2 + z^2)^0.5 - r_earth;

function [rho] = atmosphere(h)
% Standard Sea Level Values
rhoSL = 1.225;       % kg/m^3
pSL   = 101325;      % Pa
tSL   = 288.15;      % K
R     = 287.058;     % J/(kg*K)
gamma = 1.4;         

% Variables
% theta = Temperature Ratio
% sigma = Pressure Ratio
% T = Temperature
% p = Pressure
% rho = Density
% a = Speed of Sound

if h > 0 && h <= 36089 
   theta = (1-6.875*10^(-6)*h);
   T = theta*tSL;
   sigma = (1-6.875*10^(-6)*h)^5.2561;
   p = sigma*pSL;
   rho = p/(R*T);

elseif h > 36089 && h <= 65617 
   theta = 0.75189;
   T = theta*tSL;
   sigma = 0.2234*exp((4.806*10^(-5)*(36089-h)));
   p = sigma*pSL;
   rho = p/(R*T);
 
elseif h > 65617 && h <= 104990 
   theta = 0.75189+1.0577*10^(-6)*(h-65617);
   T = theta*tSL;
   sigma = 3.174716*10^(-6)*(0.75189+(1.0577*10^-6)*(h-65617))^-34.164;
   p = sigma*pSL;
   rho = p/(R*T);
end
end

%% Force Control
% Thrust Force 
F_t = 74400000;     % Constant Thrust Magnitude (Approximate Value)

% Drag Force Components 

% Velocity Squared Components 
vx = (xdot*cos(phi)*cos(theta) + ydot*sin(phi)*cos(theta) + zdot*sin(theta) + L*thetadot)^2;
vy = (ydot*cos(phi)-psidot*L*sin(theta))^2;
vz = (xdot*sin(phi)*sin(theta)+zdot*cos(theta))^2;

% Unit Vector Components 
v_mag = (vx^2 + vy^2 + vz^2)^0.5;

v1 = vx / v_mag;
v2 = vy / v_mag;
v3 = vz / v_mag;

%% Equations of Motion

% Mass Matrix 
M = [m_r*cos(phi)*cos(theta)   m_r*sin(phi)*cos(theta)   m_r*sin(theta)   0                  m_r*L;
     0                         m_r*cos(phi)              0               -m_r*L*sin(theta)   0;
     m_r*sin(phi)*sin(theta)   0                         m_r*cos(theta)  0                   0;
     0                         0                         0               Ixx*sin(theta)      0;
     0                         0                         0               0                   Iyy];

% Force Matrix  
F = [-m_r*g*sin(theta) - F_t*sin(lambda)*cos(delta) - D_mag*v1 - m_r*L*psidot^2*sin(theta)*cos(theta);
                                   - F_t*sin(lambda)*sin(delta) - D_mag*v2;
     -m_r*g*cos(theta) + F_t*cos(lambda) - D_mag*v3 + m_r*L*thetadot^2 + m_r*L*psidot^2*(sin(theta))^2;
            -L*F_t*sin(lambda)*sin(delta) - (L-c_p)*D_mag*v2 - Izz*thetadot*psidot*cos(theta);
        L*F_t*sin(lambda)*cos(delta) + (L-c_p)*D_mag*v1 -(Ixx-Izz)*psidot^2*sin(theta)*cos(theta)];

% Acceleration State Vector 2nd Order ODE
qdotdot = inv(M) * F;
