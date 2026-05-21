function [t_hist_dim, x_vec_hist_dim, Phi_mtx_hist_dim, i_e_dim, t_e_dim, x_e_vec_dim, Phi_mtx_e_dim, params_dim] ...
    = nd2dim_cr3bp(t_hist_nd, x_vec_hist_nd, Phi_mtx_hist_nd, i_e_nd, t_e_nd, x_e_vec_nd, Phi_mtx_e_nd, params_nd)
%==========================================================================
%
% Converts nondimensional trajectory data from the Circular Restricted
% Three-Body Problem (CR3BP) into dimensional units using the prescribed
% characteristic scaling.
%
% MODEL DESCRIPTION:
% The CR3BP is formulated in nondimensional units such that:
%   - Total mass: m1 + m2 = 1
%   - Distance between primaries: 1
%   - Mean motion: n = 1
%
% While this normalization simplifies the equations of motion and numerical
% integration, physical interpretation (e.g., km, km/s, seconds) requires
% mapping the nondimensional solution back to dimensional space.
%
% This function performs that mapping for:
%   - Time histories
%   - State histories (position and velocity)
%   - State Transition Matrix (STM) histories
%   - Event data (times, states, STMs)
%
% STATE DEFINITION:
% The CR3BP state is assumed to be expressed in Lagrangian coordinates:
%
%   x_vec = [x; y; z; v_x; v_y; v_z]
%
% where:
%   - Position components are nondimensionalized by the characteristic length
%   - Velocity components are nondimensionalized by L/T
%
% SCALING:
% The dimensional quantities are recovered via:
%
%   r_dim = R * r_nd
%   v_dim = (R/T) * v_nd
%   t_dim = T * t_nd
%
% where:
%   - R = characteristic length (e.g., primary separation)
%   - T = characteristic time (e.g., inverse mean motion)
%
% STM TRANSFORMATION:
% The State Transition Matrix transforms under the change of variables as:
%
%   Phi_dim = S * Phi_nd * S^{-1}
%
% where:
%
%   S = [ R*I     0
%         0     (R/T)*I ]
%
% This ensures that the STM remains consistent with the dimensionalized
% state vector.
%
% NOTES:
% - Event indices are unchanged by the transformation.
% - The structure of the trajectory is preserved; only units are changed.
% - Eigenvalues of the STM (e.g., for stability analysis) are invariant
%   under this similarity transformation.
%
%
% Author: G. Montseny
% Date: May 5, 2026
%
% INPUT:                 Description                                   Units
%
%  t_hist_nd     -   nondimensional time history (Nx1)                [-]
%  x_vec_hist_nd -   nondimensional state history (Nx6)               [-]
%  Phi_mtx_hist_nd - nondimensional STM history (Nx6x6)               [-]
%  i_e_nd        -   event indices                                    [-]
%  t_e_nd        -   nondimensional event times                       [-]
%  x_e_vec_nd    -   nondimensional event states (Ne x 6)             [-]
%  Phi_mtx_e_nd  -   nondimensional event STMs (Ne x 6 x 6)           [-]
%  params_nd     -   parameter struct (nondimensional)                [-]
%                   params_nd.scale.length  - characteristic length   [L]
%                   params_nd.scale.time    - characteristic time     [T]
%
% OUTPUT:                Description                                   Units
%
%  t_hist_dim     -   dimensional time history (Nx1)                  [T]
%  x_vec_hist_dim -   dimensional state history (Nx6)                 [L, L/T]
%  Phi_mtx_hist_dim - dimensional STM history (Nx6x6)                 [-]
%  i_e_dim        -   event indices                                   [-]
%  t_e_dim        -   dimensional event times                         [T]
%  x_e_vec_dim    -   dimensional event states (Ne x 6)               [L, L/T]
%  Phi_mtx_e_dim  -   dimensional event STMs (Ne x 6 x 6)             [-]
%  params_dim     -   parameter struct with updated units flag        [-]
%                     params_dim.model.units = 'dim'
%
%==========================================================================

    % Initialize
    t_hist_dim = zeros(size(t_hist_nd));
    x_vec_hist_dim = zeros(size(x_vec_hist_nd));
    Phi_mtx_hist_dim = zeros(size(Phi_mtx_hist_nd));
    i_e_dim = zeros(size(i_e_nd));
    t_e_dim = zeros(size(t_e_nd));
    x_e_vec_dim = zeros(size(x_e_vec_nd));
    Phi_mtx_e_dim = zeros(size(Phi_mtx_e_nd));
    N = length(t_hist_nd);
    Ne = length(t_e_nd);

    % Variables that do not change with the dimensionalization
    i_e_dim = i_e_nd;

    % Recovering scale factors
    T = params_nd.scale.time;
    R = params_nd.scale.length;
    S_mtx = [R*eye(3,3), 0*eye(3,3);
        0*eye(3,3), (R/T)*eye(3,3);];
    S_inv_mtx = [(1/R)*eye(3,3), 0*eye(3,3);
    0*eye(3,3), (T/R)*eye(3,3);];

    % Dimensionalizing time 
    t_hist_dim = T*t_hist_nd;
    t_e_dim = T*t_e_nd;

    % Dimensionalize distances and velocities
    for i = 1 : N
        x_vec_dim = S_mtx*x_vec_hist_nd(i, :)';
        x_vec_hist_dim(i, :) = x_vec_dim';
    end

    for i = 1 : Ne
        x_e_vec = S_mtx*x_e_vec_nd(i, :)';
        x_e_vec_dim(i, :) = x_e_vec';
    end

    % Dimentoinalizing the STM
    for i = 1 : N
        Phi_mtx_nd = squeeze(Phi_mtx_hist_nd(i, :,:));
        Phi_mtx_hist_dim(i, :,:) = S_mtx*Phi_mtx_nd*S_inv_mtx;
    end

    for i = 1 : Ne
        Phi_mtx_nd = squeeze(Phi_mtx_e_nd(i, :,:));
        Phi_mtx_e_dim(i, :,:) = S_mtx*Phi_mtx_nd*S_inv_mtx;
    end

    % Changing the model of the parameters
    params_dim = params_nd;
    params_dim.model.units = 'dim';

end