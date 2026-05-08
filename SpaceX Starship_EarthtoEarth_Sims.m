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
L = 55.74;                       % Distance from Base of Rocket to Center of Mass (m)
c_p = 55.09;                     % Location of Center of Pressure (m)
CD = 0.6;                        % Coefficient of Drag
A = 63.61;                       % Frontal Nose Area for Drag (m^2)
D_const = (0.5) * CD * rho * A;  % Drag Coefficient 

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
Icone = [Ixx_Co      0         0
         0         Iyy_Co      0
         0           0       Izz_Co];

Icylinder = [Ixx_Cy      0         0
             0         Iyy_Cy      0
             0           0       Izz_Cy];

% Inertia Dyadic of Rocket 
I = Icone + Icylinder;
Ixx = I(1,1);
Iyy = I(2,2);
Izz = I(3,3);

%% Atmospheric Model 



%% Drag Froce 


%% Equations of Motion
