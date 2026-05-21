function E_hist = energy_cr3bp(t_hist, x_vec_hist, params)
%==========================================================================
%
% Computes the energy-like integral of motion for the dimensionless Circular
% Restricted Three-Body Problem (CR3BP) in the uniformly rotating synodic
% frame, using the Lagrangian (position–velocity) formulation.
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
% STATE DEFINITION (LAGRANGIAN FORM):
%   x_vec = [x; y; z; v_x; v_y; v_z]
%
% where (x,y,z) is position and (v_x,v_y,v_z) is velocity in the synodic frame.
%
% ENERGY DEFINITION:
% The energy-like quantity is defined as:
%
%   E = T - V
%
% where:
%   T = 1/2*(v_x^2 + v_y^2 + v_z^2)
%   V = 1/2*(x^2 + y^2) + (1-mu)/r_13 + mu/r_23
%
% and r_13, r_23 are the distances to the primaries.
%
% This quantity is related to the Jacobi constant C by:
%
%   C = -2E
%
% NOTES:
% - The energy E is conserved along trajectories of the CR3BP.
% - This provides a useful diagnostic for numerical integration accuracy.
%
% MODEL ID: CR3BP_LAG_SYN_ND
%
% Author: G. Montseny
% Date: May 4, 2026
%
% INPUT:               Description                                   Units
%
%  t_hist      -   time vector (Nx1)                                 [-]
%  x_vec_hist  -   state history (Nx6)                               [-]
%  params      -   struct containing:                                [-]
%                    mu  - mass ratio m2/(m1+m2)                     [-]
%
% OUTPUT:              Description                                   Units
%
%  E_hist      -   energy history (Nx1)                              [-]
%
%==========================================================================

    % Initialization
    N = length(t_hist);
    E_hist = zeros(N,1);

    % Loop
    for i = 1 : N

        % Extract variables
        x_vec = x_vec_hist(i,:);
        mu = params.mu;
        x = x_vec(1); y = x_vec(2);  z = x_vec(3); 
        v_x = x_vec(4); v_y = x_vec(5); v_z = x_vec(6); 

        % Calculate kinetic and potential energy
        T = 0.5*(v_x^2+v_y^2+v_z^2);
        V = potential_energy(x, y, z);

        % Calculate total energy and append to list
        E_hist(i) = T - V;

    end

    % Helper function
    function V = potential_energy(x,y,z)
        V = 0.5*(x^2+y^2) + (1-mu)/sqrt((x+mu)^2+y^2+z^2) + mu/sqrt((x-1+mu)^2+y^2+z^2);
    end

end