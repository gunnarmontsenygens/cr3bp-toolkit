function params = params_EM_cr3bp()

%==========================================================================
%
% Builds parameter struct for the Earth–Moon Circular Restricted Three-Body
% Problem (CR3BP) in the synodic frame, using nondimensional units.
%
% MODEL ID: CR3BP_LAG_SYN_ND
%
% Author: G. Montseny
% Date: May 4, 2026
%
% OUTPUT:              Description                                   Units
%
%  params     -   parameter struct                                  [-]
%
%==========================================================================

    %-------------------------------
    % Physical constants (SI)
    %-------------------------------
    m_E = 5.97219e24;      % Earth mass          [kg]
    m_M = 7.34767309e22;   % Moon mass           [kg]

    r_EM = 384400e3;       % Earth–Moon distance [m]

    %-------------------------------
    % CR3BP mass parameter
    %-------------------------------
    mu = m_M / (m_E + m_M);

    %-------------------------------
    % Model definition
    %-------------------------------
    params.model.name = 'CR3BP';
    params.model.formulation = 'lagrangian';
    params.model.frame = 'synodic';
    params.model.units = 'nd';

    %-------------------------------
    % Core parameters
    %-------------------------------
    params.mu = mu;

    %-------------------------------
    % Functions
    %-------------------------------
    params.fun.eom = @eom_cr3bp;
    params.fun.integrate = @integrate_cr3bp;

    %-------------------------------
    % Scaling (useful later)
    %-------------------------------
    params.scale.length = r_EM;          % [m]
    params.scale.mass = m_E + m_M;       % [kg]

    % Mean motion (dimensional)
    G = 6.67430e-11;
    params.scale.n = sqrt(G*(m_E+m_M)/r_EM^3); % [rad/s]

    % Time scale
    params.scale.time = 1 / params.scale.n;    % [s]

    %-------------------------------
    % ODE options
    %-------------------------------
    params.ode.options = odeset( ...
        'RelTol', 1e-12, ...
        'AbsTol', 1e-12);

end