function dX_dt_vec = eom_ext_cr3bp(t, X_vec, params)
%==========================================================================
%
% Computes the dimensionless Circular Restricted Three-Body Problem (CR3BP)
% equations of motion in the uniformly rotating synodic frame together with
% the State Transition Matrix (STM) dynamics, using the Lagrangian
% (position–velocity) formulation.
%
% MODEL DESCRIPTION:
% The CR3BP describes the motion of a massless particle under the
% gravitational influence of two primary bodies (m1, m2) that orbit their
% barycenter in circular orbits. In the synodic frame, the primaries are
% fixed on the x-axis.
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
% The physical state evolves according to:
%
%   x_dot = f(x)
%
% and the STM evolves according to:
%
%   Phi_dot = A(x) * Phi
%
% where:
%   - f(x) is the CR3BP vector field
%   - A(x) = df/dx is the Jacobian matrix
%
% NOTES:
% - This is the LAGRANGIAN formulation (not Hamiltonian).
% - The system is autonomous; A depends on state but not explicitly on time.
% - The STM captures the linearized flow and is used for:
%     • stability analysis
%     • monodromy matrix computation
%     • invariant manifolds
%
% MODEL ID: CR3BP_LAG_SYN_ND
%
% Author: G. Montseny
% Date: May 4, 2026
%
% INPUT:               Description                                   Units
%
%  t         -   time (unused, included for ODE solver compatibility) [-]
%  X_vec     -   augmented state (42x1)                              [-]
%  params    -   struct containing:                                  [-]
%                   mu  - mass ratio m2/(m1+m2)                      [-]
%
% OUTPUT:              Description                                   Units
%
%  dX_dt_vec -   time derivative of augmented state (42x1)           [-]
%
%==========================================================================

    % Initialization
    X_vec = X_vec(:);
    
    % Extract state and STM vector
    x_vec = X_vec(1:6);
    Phi_vec = X_vec(7:42);

    % EoM
    dx_dt_vec = eom_cr3bp(t, x_vec, params);

    % STM
    Phi_mtx = reshape(Phi_vec, 6, 6);
    A_t = jacobian_cr3bp(t, x_vec, params);
    dPhi_dt_mtx = A_t*Phi_mtx;

    % Put vectors back into Y
    dx_dt_vec = dx_dt_vec(:);
    dPhi_dt_vec = dPhi_dt_mtx(:);
    dX_dt_vec = [dx_dt_vec; dPhi_dt_vec];
end