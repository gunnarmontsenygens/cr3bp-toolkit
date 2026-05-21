function A = jacobian_cr3bp(t, x_vec, params)
%==========================================================================
%
% Computes the Jacobian matrix for the dimensionless Circular Restricted
% Three-Body Problem (CR3BP) equations of motion in the uniformly rotating
% synodic frame, using the Lagrangian (position-velocity) formulation.
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
% STATE DEFINITION (LAGRANGIAN FORM):
%   x_vec = [x; y; z; v_x; v_y; v_z]
%
% where (x,y,z) is position and (v_x,v_y,v_z) is velocity in the synodic
% frame.
%
% JACOBIAN DEFINITION:
% The Jacobian matrix is defined as:
%
%   A = df/dx
%
% where f is the first-order CR3BP vector field:
%
%   d/dt(x_vec) = f(x_vec)
%
% This matrix is used to propagate the State Transition Matrix (STM):
%
%   Phi_dot = A*Phi
%
% The lower-left 3x3 block of A is the Hessian of the effective potential:
%
%   V(x,y,z) = 1/2*(x^2 + y^2)
%              + (1-mu)/r_13
%              + mu/r_23
%
% NOTES:
% - This is the Lagrangian, not Hamiltonian, CR3BP Jacobian.
% - The system is autonomous in the synodic frame, so A depends on the
%   current state but not explicitly on time.
% - Coriolis terms appear as constant velocity partial derivatives.
%
% MODEL ID: CR3BP_LAG_SYN_ND
%
% Author: G. Montseny
% Date: May 4, 2026
%
% INPUT:               Description                                   Units
%
%  t         -   time (unused, included for ODE solver compatibility) [-]
%  x_vec     -   state vector [x;y;z;v_x;v_y;v_z] (6x1)               [-]
%  params    -   struct containing:                                   [-]
%                   mu  - mass ratio m2/(m1+m2)                       [-]
%
% OUTPUT:              Description                                   Units
%
%  A         -   Jacobian matrix df/dx (6x6)                          [-]
%
%==========================================================================
    
    % Initialization
    x_vec = x_vec(:);

    % Extract params
    mu = params.mu;
    x_1 = -mu; r_1_vec = [x_1, 0, 0]';
    x_2 = 1-mu; r_2_vec = [x_2, 0, 0]';

    % Extract values from x_vec
    x = x_vec(1);
    y = x_vec(2);
    z = x_vec(3);
    
    % Calculate important quantities
    r_vec = [x; y; z];
    r_13 = norm(r_vec-r_1_vec);
    r_23 = norm(r_vec-r_2_vec);

    % Potential V = 0.5*(x^2+y^2) + U(x,y,z) double partial derivatives
    V_xx = 1 + 3*(1-mu)*(x-x_1)^2/r_13^5 + 3*mu*(x-x_2)^2/r_23^5 - (1-mu)/r_13^3 - mu/r_23^3;
    V_xy = 3*(1-mu)*(x-x_1)*y/r_13^5 + 3*mu*(x-x_2)*y/r_23^5;
    V_xz = 3*(1-mu)*(x-x_1)*z/r_13^5 + 3*mu*(x-x_2)*z/r_23^5;
    V_yy = 1 + 3*(1-mu)*y^2/r_13^5 + 3*mu*y^2/r_23^5 - (1-mu)/r_13^3 - mu/r_23^3;
    V_yz = 3*(1-mu)*y*z/r_13^5 + 3*mu*y*z/r_23^5;
    V_zz = 3*(1-mu)*z^2/r_13^5 + 3*mu*z^2/r_23^5 - (1-mu)/r_13^3 - mu/r_23^3;


    % Dynamics matrix A
    A = [0,    0,    0,    1,  0, 0;
         0,    0,    0,    0,  1, 0;
         0,    0,    0,    0,  0, 1;
         V_xx, V_xy, V_xz, 0,  2, 0;
         V_xy, V_yy, V_yz, -2, 0, 0;
         V_xz, V_yz, V_zz, 0,  0, 0];
        
end










