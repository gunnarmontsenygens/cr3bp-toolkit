function [t_hist, x_vec_hist, Phi_mtx_hist, i_e, t_e, x_e_vec, Phi_mtx_e] = integrate_cr3bp(t_span, x_0_vec, params, event_fun)
%==========================================================================
%
% Integrates the dimensionless Circular Restricted Three-Body Problem (CR3BP)
% equations of motion in the uniformly rotating synodic frame together with
% the State Transition Matrix (STM) dynamics, using the Lagrangian
% (position–velocity) formulation.
%
% MODEL DESCRIPTION:
% The CR3BP describes the motion of a massless particle under the
% gravitational influence of two primary bodies (m1, m2) that orbit their
% barycenter in circular orbits. In the synodic frame, the primaries are
% fixed on the x-axis, resulting in an autonomous dynamical system.
%
% NORMALIZATION:
% The system is nondimensionalized such that:
%   - Total mass: m1 + m2 = 1
%   - Distance between primaries: 1
%   - Mean motion: n = 1
%
% In this frame:
%   Primary 1 is located at x = -mu
%   Primary 2 is located at x = 1 - mu
%
% STATE DEFINITION (AUGMENTED SYSTEM):
% The augmented state combines the physical state and the STM:
%
%   X_vec = [x; y; z; v_x; v_y; v_z; vec(Phi)]
%
% where:
%   - x_vec = [x; y; z; v_x; v_y; v_z] is the 6x1 Lagrangian state
%   - Phi is the 6x6 State Transition Matrix
%   - vec(Phi) is the column-wise vectorization of Phi (36x1)
%
% Total state dimension: 6 + 36 = 42
%
% DYNAMICS:
% The system is integrated by solving:
%
%   x_dot   = f(x)
%   Phi_dot = A(x) * Phi
%
% where:
%   - f(x) is the CR3BP vector field
%   - A(x) = df/dx is the Jacobian matrix
%
% The STM is initialized as:
%
%   Phi(0) = I
%
% EVENT HANDLING:
% An optional event function can be provided to detect specific conditions
% during integration (e.g., plane crossings, impacts, section conditions).
% If provided, the integrator returns:
%   - event times
%   - event states
%   - corresponding STM values
%
% NOTES:
% - This is the LAGRANGIAN formulation (not Hamiltonian).
% - The system is autonomous; dynamics depend on state but not explicitly on time.
% - The STM captures the linearized flow and is used for:
%     • stability analysis (eigenvalues of monodromy matrix)
%     • invariant manifold computation
%     • sensitivity analysis and targeting
%
% MODEL ID: CR3BP_LAG_SYN_ND
%
% Author: G. Montseny
% Date: May 4, 2026
%
% INPUT:                 Description                                   Units
%
%  t_span      -   time span for integration [t0 tf]                   [-]
%  x_0_vec     -   initial state vector (6x1)                          [-]
%  params      -   parameter struct                                    [-]
%                  params.mu          - mass ratio m2/(m1+m2)          [-]
%                  params.ode.options - ODE solver options             [-]
%  event_fun   -   optional event function handle                      [-]
%
% OUTPUT:                Description                                   Units
%
%  t_hist       -   time vector (Nx1)                                  [-]
%  x_vec_hist   -   state history (Nx6)                                [-]
%  Phi_mtx_hist -   STM history (Nx6x6)                                [-]
%  i_e          -   event indices                                      [-]
%  t_e          -   event times                                        [-]
%  x_e_vec      -   state at events (Ne x 6)                           [-]
%  Phi_mtx_e    -   STM at events (Ne x 6 x 6)                         [-]
%
%==========================================================================


    % Check if there's any events
    if nargin < 4
        event_fun = [];
    end
    
    % Extract base ODE options from params 
    ode_options = params.ode.options;
    
    % Implement event in case it's not empy
    if ~isempty(event_fun)
        ode_options = odeset(ode_options, 'Events', @(t,Y) event_fun(t,Y,params));
    end
    
    % Initialization
    x_0_vec = x_0_vec(:);
    n = length(x_0_vec);
    Phi_0_mtx = eye(n);
    X_0_vec = [x_0_vec; Phi_0_mtx(:)];
    
    % Integration
    if isempty(event_fun)

        [t_hist, X_vec_hist] = ode113( @(t,X) eom_ext_cr3bp(t, X, params), t_span, X_0_vec, ode_options);

        % Empty event outputs
        t_e = [];
        x_e_vec = [];
        Phi_mtx_e = [];
        i_e = [];

    else
        [t_hist, X_vec_hist, t_e, X_e_vec, i_e] = ode113(@(t,X) eom_ext_cr3bp(t, X, params), t_span, X_0_vec, ode_options);
    
        x_e_vec = X_e_vec(:,1:n);
        N_e = size(X_e_vec,1);
        Phi_mtx_e = zeros(N_e,n,n);
    
        for k = 1:N_e
            Phi_mtx_e(k,:,:) = reshape(X_e_vec(k,n+1:n+n^2), n, n);
        end
    end

    % Recover state history
    x_vec_hist = X_vec_hist(:,1:n);

    % Recover STM history as (N x n x n)
    N = size(X_vec_hist,1);
    Phi_mtx_hist = zeros(N,n,n);
    
    for k = 1:N
        Phi_mtx_hist(k,:,:) = reshape(X_vec_hist(k,n+1:n + n^2),n,n);
    end
end