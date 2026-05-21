function dx_dt_vec =  eom_cr3bp(t, x_vec, params)
%==========================================================================
%
% Computes the dimensionless Circular Restricted Three-Body Problem (CR3BP)
% equations of motion in the uniformly rotating synodic frame, using the
% Lagrangian (position–velocity) formulation.
%
% MODEL DESCRIPTION:
% The CR3BP describes the motion of a massless particle under the
% gravitational influence of two primary bodies (m1, m2) that orbit their
% barycenter in circular orbits. The equations are written in the synodic
% (rotating) frame where the primaries are fixed on the x-axis.
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
% EQUATIONS OF MOTION:
%
%   x_ddot =  2*y_dot + dV/dx
%   y_ddot = -2*x_dot + dV/dy
%   z_ddot =           dV/dz
%
% where the effective potential is:
%
%   V(x,y,z) = 1/2*(x^2 + y^2)
%              + (1-mu)/r_13
%              + mu/r_23
%
% and r_13, r_23 are the distances to the primaries.
%
% NOTES:
% - This is the LAGRANGIAN (position–velocity) formulation of the CR3BP.
% - The system is autonomous in the synodic frame.
% - Coriolis and centrifugal terms arise from the rotating frame.
% - The equivalent Hamiltonian formulation can be obtained via:
%       p_x = v_x - y
%       p_y = v_y + x
%
% MODEL ID: CR3BP_LAG_SYN_ND
%
% Author: G. Montseny
% Date: May 4, 2026
%
% INPUT:               Description                                   Units
%
%  t         -   time (unused, included for ODE solver compatibility) [-]
%  x_vec     -   state vector (6x1)                                  [-]
%  params    -   struct containing:                                  [-]
%                   mu  - mass ratio m2/(m1+m2)                      [-]
%
% OUTPUT:              Description                                   Units
%
%  dx_dt_vec -   time derivative of state (6x1)                      [-]
%
%==========================================================================

    % Initialization
    x_vec = x_vec(:);
    dx_dt_vec = zeros(size(x_vec));

    % Extract params
    mu = params.mu;
    x_1 = -mu; r_1_vec = [x_1, 0, 0]';
    x_2 = 1-mu; r_2_vec = [x_2, 0, 0]';

    % Extract values from x_vec
    x = x_vec(1);
    y = x_vec(2);
    z = x_vec(3);
    v_x = x_vec(4);
    v_y = x_vec(5);
    v_z = x_vec(6);
    
    % Calculate important quantities
    r_vec = [x; y; z];
    r_13 = norm(r_vec-r_1_vec);
    r_23 = norm(r_vec-r_2_vec);

    % Potential V = 0.5*(x^2+y^2) + U(x,y,z) derivatives
    dVdx = x - (1-mu)*(x-x_1)/r_13^3 - mu*(x-x_2)/r_23^3;
    dVdy = y - (1-mu)*y/r_13^3 - mu*y/r_23^3;
    dVdz = - (1-mu)*z/r_13^3 - mu*z/r_23^3;

    % Lagrangian EoM
    dx_dt_vec(1) = v_x;
    dx_dt_vec(2) = v_y;
    dx_dt_vec(3) = v_z;
    dx_dt_vec(4) =  2*v_y + dVdx;
    dx_dt_vec(5) = -2*v_x + dVdy;
    dx_dt_vec(6) = dVdz;

end