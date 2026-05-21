function L_pts = lagrange_points_cr3bp(mu)
%==========================================================================
%
% Computes the Lagrange (equilibrium) points of the dimensionless Circular
% Restricted Three-Body Problem (CR3BP) in the uniformly rotating synodic
% frame.
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
% LAGRANGE POINTS:
% The five equilibrium points are defined as the solutions to:
%
%   dV/dx = 0,  dV/dy = 0,  dV/dz = 0
%
% where V(x,y,z) is the effective potential of the CR3BP.
%
% - The collinear points (L1, L2, L3) lie along the x-axis and are computed
%   numerically by solving dV/dx = 0 with y = z = 0.
%
% - The triangular points (L4, L5) form equilateral triangles with the
%   primaries and have analytical positions:
%
%       L4 = [1/2 - mu,  sqrt(3)/2]
%       L5 = [1/2 - mu, -sqrt(3)/2]
%
% NOTES:
% - L1, L2, L3 are saddle-type equilibrium points.
% - L4 and L5 are linearly stable for mu < mu_crit ≈ 0.0385.
% - These points are widely used as reference locations for periodic orbits
%   and invariant manifold computations.
%
% MODEL ID: CR3BP_LAG_SYN_ND
%
% Author: G. Montseny
% Date: May 4, 2026
%
% INPUT:               Description                                   Units
%
%  mu        -   mass ratio m2/(m1+m2)                               [-]
%
% OUTPUT:              Description                                   Units
%
%  L_pts     -   struct containing Lagrange points                   [-]
%                L1, L2, L3, L4, L5 as 6D positions 
%
%==========================================================================

    % Initialization
    options = optimset('TolX',1e-12,'TolFun',1e-12);

    % Collinear points
    eps = 1e-6;
    x_L1 = fzero(@(x) V_x(x,0,0), [0.5 - mu, 1 - mu - eps], options);
    x_L2 = fzero(@(x) V_x(x,0,0), [1 - mu + eps, 2], options);
    x_L3 = fzero(@(x) V_x(x,0,0), [-2, - mu - eps], options);

    % Triangular points
    L4 = [0.5 - mu,  sqrt(3)/2, 0, 0, 0, 0];
    L5 = [0.5 - mu, -sqrt(3)/2, 0, 0, 0, 0];

    % Store output
    L_pts.L1 = [x_L1, 0, 0, 0, 0, 0];
    L_pts.L2 = [x_L2, 0, 0, 0, 0, 0];
    L_pts.L3 = [x_L3, 0, 0, 0, 0, 0];
    L_pts.L4 = L4;
    L_pts.L5 = L5;

    function dVdx = V_x(x,y,z)
        dVdx = x - (1-mu)*(x+mu)/((x+mu)^2 + y^2 + z^2)^(3/2) - mu*(x-1+mu)/((x-1+mu)^2 + y^2 + z^2)^(3/2);
    end

end