function dCjdx_vec = jacobi_gradient_cr3bp(x_vec, params)
%==========================================================================
%
% Computes the gradient of the Jacobi constant for the dimensionless
% Circular Restricted Three-Body Problem (CR3BP) in the uniformly rotating
% synodic frame, using the Lagrangian (position–velocity) formulation.
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
% JACOBI CONSTANT:
% The Jacobi constant is defined as:
%
%   C_J = 2*Omega(x,y,z) - (v_x^2 + v_y^2 + v_z^2)
%
% where:
%
%   Omega(x,y,z) = 1/2*(x^2 + y^2)
%                  + (1-mu)/r_1
%                  + mu/r_2
%
% and:
%
%   r_1 = sqrt((x+mu)^2 + y^2 + z^2)
%   r_2 = sqrt((x-1+mu)^2 + y^2 + z^2)
%
% The gradient of the Jacobi constant is:
%
%   grad(C_J) = [2*Omega_x;
%                2*Omega_y;
%                2*Omega_z;
%               -2*v_x;
%               -2*v_y;
%               -2*v_z]
%
% NOTES:
% - The Jacobi gradient is commonly used in differential correction and
%   continuation methods for periodic orbit computation.
% - This implementation assumes the Lagrangian CR3BP state convention.
%
% MODEL ID: CR3BP_LAG_SYN_ND
%
% Author: G. Montseny
% Date: May 6, 2026
%
% INPUT:               Description                                   Units
%
%  x_vec      -   state vector [x;y;z;v_x;v_y;v_z]                  [-]
%  params     -   struct containing:                                [-]
%                   mu  - mass ratio m2/(m1+m2)                     [-]
%
% OUTPUT:              Description                                   Units
%
%  dCjdx_vec  -   gradient of Jacobi constant wrt x_vec             [-]
%
%==========================================================================

    % Extract parameters
    mu = params.mu;

    % Extract state variables
    x   = x_vec(1);
    y   = x_vec(2);
    z   = x_vec(3);
    v_x = x_vec(4);
    v_y = x_vec(5);
    v_z = x_vec(6);

    % Distances to primaries
    r_1 = sqrt((x + mu)^2 + y^2 + z^2);
    r_2 = sqrt((x - 1 + mu)^2 + y^2 + z^2);

    % Derivatives of pseudo-potential Omega
    Omega_x = x ...
        - (1-mu)*(x + mu)/r_1^3 ...
        - mu*(x - 1 + mu)/r_2^3;

    Omega_y = y ...
        - (1-mu)*y/r_1^3 ...
        - mu*y/r_2^3;

    Omega_z = ...
        - (1-mu)*z/r_1^3 ...
        - mu*z/r_2^3;

    % Gradient of Jacobi constant
    dCjdx_vec = [ 2*Omega_x;
                  2*Omega_y;
                  2*Omega_z;
                 -2*v_x;
                 -2*v_y;
                 -2*v_z ];

end