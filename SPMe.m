%% SPMe Model Implementation in MATLAB
syms t

%% ===========================
% 1. Global Parameters
% ===========================
global N M NM f y0 Numexp Totexp

% C-rate
Crate = -1;   % negative = discharge

% Experimental data
Numexp = 63;        % Number of experimental data points
Totexp = 3100;      % Total discharge time (s)

%% ===========================
% 2. Discretization
% ===========================
% Node number for electrolyte (positive, separator, negative)
N  = 2;   % positive electrode
M  = 2;   % separator
NM = 2;   % negative electrode

%% ===========================
% 3. Design Parameters
% ===========================
ep    = 0.335;    % Porosity at positive
es    = 0.47;     % Porosity at separator
en    = 0.25;     % Porosity at negative
brugp = 2.43;
brugs = 2.57;
brugn = 2.91;

lp  = 75.6e-6;    % Thickness positive electrode
ls  = 12e-6;      % Thickness separator
ln  = 85.2e-6;    % Thickness negative electrode

Rpp = 5.22e-6;    % Particle radius positive
Rpn = 5.86e-6;    % Particle radius negative

F = 96487;        % Faraday constant (C/mol)
R = 8.3143;       % Gas constant (J/(mol K))
T = 298.15;       % Temperature (K)
tplus = 0.363;    % Transference number

ap = (3/Rpp)*(1-ep);   % Specific interfacial area (pos)
an = (3/Rpn)*(1-en);   % Specific interfacial area (neg)

Acell = 0.11;       % Electrode area (m^2)
Capa  = 5;          % Nominal capacity (Ah)
iapp  = Capa * Crate / Acell; % Applied current density (A/m^2)

%% ===========================
% 4. Transport and Kinetic Parameters
% ===========================
c0     = 1000;         % Initial electrolyte concentration (mol/m^3)
D1     = 1e-9;         % Electrolyte diffusion coefficient (m^2/s)
Kappa  = 1.17;         % Electrolyte conductivity (S/m)
ctp    = 51765;        % Max solid concentration positive
ctn    = 29583;        % Max solid concentration negative
sigmap = 0.18;         % Solid conductivity positive (S/m)
sigman = 215;          % Solid conductivity negative (S/m)
Dsp    = 4e-15;        % Solid diffusion positive (m^2/s)
Dsn    = 3.3e-14;      % Solid diffusion negative (m^2/s)
kp     = 0.7e-11;      % Reaction rate positive
kn     = 0.7e-12;      % Reaction rate negative

% Effective properties
Keffp = Kappa * ep^brugp;
Keffs = Kappa * es^brugs;
Keffn = Kappa * en^brugn;

D2pos = ep^brugp * D1;
D2sep = es^brugs * D1;
D2neg = en^brugn * D1;

%% ===========================
% 5. Discretization step sizes
% ===========================
h  = lp/(N+1);   % positive
h2 = ls/(M+1);   % separator
h3 = ln/(NM+1);  % negative

%% ===========================
% 6. Define symbolic variables
% ===========================
% State variables:
% u1: electrolyte concentration (1D discretized)
% csp: average solid concentration positive
% csn: average solid concentration negative

Nt = 1 + N + 1 + M + 1 + NM + 1 + 2;  % u1 nodes + 2 solid concentrations
X  = cell(Nt,1);

for ii = 1:Nt
    X{ii} = symfun(str2sym(sprintf('X_%d(t)',ii)), t);
end

varsX = [X{:}];

% Split state variables
u1 = X(1:1+N+1+M+1+NM+1); % electrolyte concentration
csp = X{end-1};           % avg solid concentration pos
csn = X{end};             % avg solid concentration neg

%% ===========================
% 7. Electrolyte dynamics (1D PDE discretized)
% ===========================
eq1 = sym(zeros(1,1+N+1+M+1+NM+1));

% Boundary conditions (Neumann)
dudxf1 = ( -u1{3} - 3*u1{1} + 4*u1{2} ) / (2*h);
bc11 = dudxf1;          % Zero flux at x=0

eq1(1) = 0 == bc11;

% Positive electrode internal nodes
for i = 2:N+1
    d2udx2 = (u1{i-1} - 2*u1{i} + u1{i+1}) / (h^2);
    eq1(i) = diff(u1{i}) == (D2pos*d2udx2 + ap*(1-tplus)*iapp/c0/ep)/ep;
end

% Separator nodes
for i = N+2+1 : N+2+M
    d2udx2 = (u1{i-1} - 2*u1{i} + u1{i+1}) / (h2^2);
    eq1(i) = diff(u1{i}) == D2sep*d2udx2 / es;
end

% Negative electrode nodes
for i = N+2+M+1+1 : N+2+M+1+NM
    d2udx2 = (u1{i-1} - 2*u1{i} + u1{i+1}) / (h3^2);
    eq1(i) = diff(u1{i}) == (D2neg*d2udx2 - an*(1-tplus)*iapp/c0/en)/en;
end

%% ===========================
% 8. Solid concentration dynamics (single particle ODE)
% ===========================
% Using simple average concentration dynamics:
% dcs/dt = -3*j/(R_s * cmax)
jpos = iapp / (F*ap);
jneg = -iapp / (F*an);

eq2 = sym(zeros(1,2));
eq2(1) = diff(csp) == -3*jpos / (Rpp * ctp);
eq2(2) = diff(csn) == -3*jneg / (Rpn * ctn);

%% ===========================
% 9. Terminal voltage computation
% ===========================
% Compute open-circuit potentials:
theta_p = csp / ctp;
Up = -0.8090*theta_p + 4.4875 - 0.0428*tanh(18.5138*(theta_p-0.5542)) ...
    - 17.7326*tanh(15.7890*(theta_p-0.3117)) + 17.5842*tanh(15.9308*(theta_p-0.3120));

theta_n = csn / ctn;
Un = 1.9793*exp(-39.3631*theta_n) + 0.2482 - 0.0909*tanh(29.8538*(theta_n-0.1234)) ...
    - 0.04478*tanh(14.9159*(theta_n-0.2769)) - 0.0205*tanh(30.4444*(theta_n-0.6103));

eta_p = asinh(iapp/(2*kp*F*((ctp - csp)*csp)^0.5));
eta_n = asinh(iapp/(2*kn*F*((ctn - csn)*csn)^0.5));

Vcell = Up - Un + eta_p - eta_n - iapp*(lp/sigmap + ln/sigman);

%% ===========================
% 10. Assemble total equations
% ===========================
eqs = [eq1, eq2];
vars = varsX;

[MM,f] = massMatrixForm(eqs, vars);
MM = odeFunction(MM, vars);
f = odeFunction(f, vars);

%% ===========================
% 11. Initial conditions
% ===========================
U0 = ones(1, length(vars));
y0 = U0;

%% ===========================
% 12. Solve system
% ===========================
opt = odeset('Mass', MM, 'MStateDependence','weak', 'RelTol',1e-5, 'AbsTol',1e-5);
tspan = [0, Totexp];
[T,Y] = ode15s(@(t,y) f(t,y), tspan, y0, opt);

%% ===========================
% 13. Compute and plot voltage
% ===========================
V = zeros(size(T));
for k = 1:length(T)
    V(k) = double(subs(Vcell, {csp,csn}, {Y(k,end-1),Y(k,end)}));
end

figure;
plot(T, V, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Voltage (V)');
title('SPMe Simulation');
