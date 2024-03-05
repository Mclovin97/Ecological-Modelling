% Parameters for Advection and Diffusion
Param.u = 1;
Param.n = 100; % Number of grids
Param.DeltaZ = 1; % Size of each grid cell
Param.Depth = Param.n/Param.DeltaZ; % Total depth
Param.Zbottom = Param.Depth; % Bottom of the grid
Param.D = 5; % Diffusive coefficient

% Parameters for Phytoplankton
Param.I0 = 350; % Incident light Intensity
Param.gmax = 1.5; % Maximal specific production rate
Param.kp = 15*10^-12; % Light attenuation of phytoplankton
Param.kw = 0.2; % Background turbidity
Param.m = 0.01 *24; % Death rate/loss rate
Param.tz = 20; % Thermocline depth
Param.H = 30; % Half-saturation constant of light-limited growth
Param.TD = 5*60*60*24; % Turbulent diffusion
Param.VV = 0.04*24; % Vertical velocity
Param.r = Param.gmax - Param.m; % Intrinsic growth rate

% Parameters for Nutrients
Param.N0 = 100; % Initial nutrient concentration
Param.rem_rate = 0.01 * 24; % Remineralization rate of dead phytoplankton to nutrients
Param.eps = 0.5; % Dimensionless
Param.HI = 45; %µmol photons/m2/s Half saturation constand of light limited growth
Param.HN =0.025; % Half saturation constant of nutrient limitited growth
Param.V = 0.042; % m/h Sinking velocity
Param.Nb = 50; % Concentration at the bottom
Param.DU = 0.3*(60*60*24); %Turbulent diffusitivy in upper mixed layer Cm^2/day
Param.DD = 50*(60*60*24); %Turbulent diffusitivity in the deep layers Cm^2/day
Param.alpha = 1*10^-9; % Nutrient content of phytoplankton

% Set up grid
Param.Z = (1/2 * Param.DeltaZ):Param.DeltaZ:(Param.Zbottom - 1/2 * Param.DeltaZ);

% Set initial conditions for Phy at time 0
Phy = zeros(Param.n, 1); % Initial condition vector for phytoplankton
middleindex = round(Param.n/2); 
Phy(middleindex) = 100000e10; % Initial peak in the middle

Nut = Param.N0 * ones(Param.n, 1); % Initial condition vector for nutrients

InitialConditions = [Phy; Nut];


% ODE solver settings
tSpan = [0, 1000]; % Time span for the simulation
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

% Solve the ODE
[T, Sol] = ode45(@(t, y) derivative(t, y, Param), tSpan, InitialConditions, options);
 Phy = Sol(1:Param.n);
 Nut = Sol(Param.n+1:2*Param.n);


% Assuming Sol is a matrix where each row represents a time step
% and the first half of the columns represent the Phy at each depth
for i = 1:length(T)
    LightIntensity(i, :) = CalcLightIntensity(Param, Sol(i, 1:Param.n), T(i))';
end


% Assuming you want to plot the profile at the last time step
finalTimeIndex = length(T);

% Extract the final Phytoplankton and Nutrient solutions
finalPhy = Sol(finalTimeIndex, 1:Param.n);
finalNut = Sol(finalTimeIndex, Param.n+1:end);

% Calculate the final Light Intensity
finalLightIntensity = CalcLightIntensity(Param, finalPhy, T);

% Calculate the Limiting Resource
finalLimitingResource = min(finalLightIntensity ./ max(finalLightIntensity), finalNut ./ max(finalNut));

figure;
plot(finalPhy, Param.Z);
set(gca, 'YDir', 'reverse');
xlabel('Phytoplankton Concentration');
ylabel('Depth (m)');
title('Phytoplankton');

% Create a figure window
figure;

% Subplot 1: Vertical profile of Phytoplankton
subplot(2, 2, 1);
plot(finalPhy, Param.Z);
set(gca, 'YDir', 'reverse');
xlabel('Phytoplankton Concentration');
ylabel('Depth (m)');
title('Phytoplankton');

% Subplot 2: Vertical profile of Nutrients
subplot(2, 2, 2);
plot(finalNut, Param.Z);
set(gca, 'YDir', 'reverse');
xlabel('Nutrient Concentration');
ylabel('Depth (m)');
title('Nutrients');

% Subplot 3: Vertical profile of Light Intensity
subplot(2, 2, 3);
plot(finalLightIntensity, Param.Z);
set(gca, 'YDir', 'reverse');
xlabel('Light Intensity');
ylabel('Depth (m)');
title('Light Intensity');

% Subplot 4: Vertical profile of Limiting Resource
subplot(2, 2, 4);
plot(finalLimitingResource, Param.Z);
set(gca, 'YDir', 'reverse');
xlabel('Limiting Resource');
ylabel('Depth (m)');
title('Limiting Resource');

% Assuming 'Sol' is a matrix where each row represents a time step,
% and columns represent Phytoplankton followed by Nutrients at each depth.

% Initialize matrices to store Phytoplankton and Nutrient data over time
PhySol = zeros(length(T), Param.n);
NutSol = zeros(length(T), Param.n);

for i = 1:length(T)
    % Extract Phy and Nut for each time step
    PhySol(i, :) = Sol(i, 1:Param.n);
    NutSol(i, :) = Sol(i, Param.n+1:end);
end

% Now, PhySol and NutSol matrices contain the Phytoplankton and Nutrient
% data respectively, for each depth (column) and time step (row).

% Plot Nutrients over time and depth
figure;
surf(T, Param.Z, NutSol');
shading interp;
xlabel('Time');
ylabel('Depth');
zlabel('Nutrients');
title('Nutrients over time and depth');
set(gca, 'YDir', 'reverse');

%Plot light intensity
figure;
imagesc(T, Param.Z, LightIntensity');
colorbar;
caxis([min(LightIntensity(:)), max(LightIntensity(:))]); % Adjust the color axis to the data range
xlabel('Time');
ylabel('Depth');
title('Light Intensity over Depth and Time');
set(gca, 'YDir', 'reverse');

function I = CalcLightIntensity (Param, Phy, t)
    %Seasonal aspect
    season = sin(2 * pi * t / 365 - pi/2); % varies between 0 and 1 over the year
    I0_baseline=Param.I0; % Average incident light intensity
    I0_amplitude=Param.I0/2; % Amplitude of seasonal variation
    I0 = I0_baseline + I0_amplitude *season;

    % Calculate I0 based on time
    % Using the trapezoidal rule for numerical integration
    Integral = cumsum((Param.kw + Param.kp * Phy) * Param.DeltaZ + ...
                      0.5 * Param.DeltaZ * (Param.kw + Param.kp * Phy));
    I = I0 * exp(-Integral);
end


function dYdt = derivative(t, Y, Param)
    % Extracting solutions for each component
    Phy = Y(1:Param.n);
    Nut = Y(Param.n+1:2*Param.n);

    % Initialize fluxes
    Ja = zeros(Param.n, 1); % Advection flux for Phy
    Jd = zeros(Param.n, 1); % Diffusion flux for Phy
    JaN = zeros(Param.n, 1); % Advection flux for Nut
    JdN = zeros(Param.n, 1); % Diffusion flux for Nut

    % Calculate advection and diffusion fluxes for Phy and Nut
    for ix = 2:Param.n
        Ja(ix) = Param.V * Phy(ix-1);
        Jd(ix) = -Param.D * (Phy(ix) - Phy(ix-1)) / Param.DeltaZ;
        JaN(ix) = Param.V * Nut(ix-1);
        JdN(ix) = -Param.D * (Nut(ix) - Nut(ix-1)) / Param.DeltaZ;
    end

    % Boundary conditions
    % Phytoplankton - no flux at the bottom, no advection at the surface
    Ja(1) = 0;
    Jd(1) = 0;
    Ja(Param.n) = 0;
    Jd(Param.n) = 0;
    
    % Nutrients - fixed concentration at the bottom, no advection at the surface
    JaN(1) = 0;
    JdN(1) = 0;
    JaN(Param.n) = 0;
    JdN(Param.n) = -Param.D .* (Param.Nb - Nut(end)) / Param.DeltaZ;

    % Combine advection and diffusion
    J = Ja + Jd;
    JN = JaN + JdN;

    % Advection-diffusion term for Phy and Nut
    advec_diffP = (-(J(3:Param.n) - J(2:Param.n-1)) / Param.DeltaZ);
    advec_diffN = (-(JN(3:Param.n) - JN(2:Param.n-1)) / Param.DeltaZ);

    % Pad advec_diffP and advec_diffN with zeros for the surface and bottom boundary
    advec_diffP = [0; advec_diffP; 0];
    advec_diffN = [0; advec_diffN; 0];

    % Light intensity at each depth (for growth rate calculation)
    I = CalcLightIntensity(Param, Phy,t);

    % Growth rate limited by light and nutrient availability
    II = I ./ (I + Param.HI);
    NN = Nut ./ (Nut + Param.HN);
    Lim = min(II, NN);
    g = Param.gmax .* Lim;

    % Change in Phytoplankton over time
    dPhy_dt = advec_diffP + g .* Phy - Param.m * Phy;

    % Change in Nutrients over time, including the remineralization
    dNut_dt = advec_diffN + (Param.m * Phy * Param.rem_rate) - (g .* Phy * Param.alpha);

    % Combine derivatives into a column vector for the ODE solver
    dYdt = [dPhy_dt; dNut_dt];
end
