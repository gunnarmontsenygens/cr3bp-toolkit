function C_hist = jacobi_constant_cr3bp(t_hist, x_vec_hist, params)
%==========================================================================
%
% Computes the Jacobi constant for the dimensionless Circular Restricted
% Three-Body Problem (CR3BP) in the uniformly rotating synodic frame, using
% the Lagrangian (position–velocity) formulation.
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
%   C = -2E
%
% where the energy-like quantity E is:
%
%   E = 1/2*(v_x^2 + v_y^2 + v_z^2) - V(x,y,z)
%
% with:
%
%   V(x,y,z) = 1/2*(x^2 + y^2)
%              + (1-mu)/r_13
%              + mu/r_23
%
% and r_13, r_23 are the distances to the primaries.
%
% NOTES:
% - The Jacobi constant is conserved along trajectories of the CR3BP.
% - It defines zero-velocity surfaces and allowed regions of motion.
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
%  C_hist      -   Jacobi constant history (Nx1)                     [-]
%
%==========================================================================

C_hist = -2*energy_cr3bp(t_hist, x_vec_hist, params);