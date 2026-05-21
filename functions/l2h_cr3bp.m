function [t_hist_h, x_vec_hist_h, Phi_mtx_hist_h, i_e_h, t_e_h, x_e_vec_h, Phi_mtx_e_h, params_h] ...
    = l2h_cr3bp(t_hist_l, x_vec_hist_l, Phi_mtx_hist_l, i_e_l, t_e_l, x_e_vec_l, Phi_mtx_e_l, params_l)
%==========================================================================
%
% Converts trajectory data for the Circular Restricted Three-Body Problem
% (CR3BP) from Lagrangian coordinates to Hamiltonian coordinates.
%
% MODEL DESCRIPTION:
% The CR3BP may be formulated using either Lagrangian position-velocity
% coordinates or Hamiltonian canonical coordinates. In the rotating synodic
% frame, these two representations are related by a linear transformation
% between velocities and canonical momenta.
%
% This function maps trajectory data from the Lagrangian state definition
%
%   x_l_vec = [x; y; z; v_x; v_y; v_z]
%
% to the Hamiltonian state definition
%
%   x_h_vec = [x; y; z; p_x; p_y; p_z]
%
% using the rotating-frame momentum relations
%
%   p_x = v_x - y
%   p_y = v_y + x
%   p_z = v_z
%
% The transformation is applied to:
%   - Time histories
%   - State histories
%   - State Transition Matrix (STM) histories
%   - Event data
%   - Parameter structure formulation flag
%
% STATE TRANSFORMATION:
% The coordinate transformation is written as
%
%   x_h_vec = T * x_l_vec
%
% where
%
%   T = [ I   0
%         S   I ]
%
% and
%
%   S = [ 0  -1   0
%         1   0   0
%         0   0   0 ]
%
% STM TRANSFORMATION:
% Since the Lagrangian-to-Hamiltonian map is linear and time-independent,
% the State Transition Matrix transforms by similarity:
%
%   Phi_h = T * Phi_l * T^{-1}
%
% This ensures that the STM remains consistent with the transformed state
% vector and maps Hamiltonian-coordinate perturbations correctly.
%
% NOTES:
% - Time histories and event times are unchanged by this transformation.
% - Event indices are unchanged.
% - Eigenvalues of the STM are invariant under this similarity
%   transformation.
% - This transformation is valid for the standard rotating-frame CR3BP
%   Lagrangian and Hamiltonian coordinate definitions.
%
%
% Author: G. Montseny
% Date: May 5, 2026
%
% INPUT:                 Description                                   Units
%
%  t_hist_l      -   time history in Lagrangian formulation (Nx1)     [-]
%  x_vec_hist_l  -   Lagrangian state history (Nx6)                   [-]
%  Phi_mtx_hist_l -  Lagrangian STM history (Nx6x6)                   [-]
%  i_e_l         -   event indices                                    [-]
%  t_e_l         -   event times                                      [-]
%  x_e_vec_l     -   Lagrangian event states (Ne x 6)                 [-]
%  Phi_mtx_e_l   -   Lagrangian event STMs (Ne x 6 x 6)               [-]
%  params_l      -   parameter struct in Lagrangian formulation       [-]
%
% OUTPUT:                Description                                   Units
%
%  t_hist_h      -   time history in Hamiltonian formulation (Nx1)    [-]
%  x_vec_hist_h  -   Hamiltonian state history (Nx6)                  [-]
%  Phi_mtx_hist_h -  Hamiltonian STM history (Nx6x6)                  [-]
%  i_e_h         -   event indices                                    [-]
%  t_e_h         -   event times                                      [-]
%  x_e_vec_h     -   Hamiltonian event states (Ne x 6)                [-]
%  Phi_mtx_e_h   -   Hamiltonian event STMs (Ne x 6 x 6)              [-]
%  params_h      -   parameter struct with updated formulation flag   [-]
%                    params_h.model.formulation = 'hamiltonian'
%
%==========================================================================

    % Initialize
    t_hist_h = zeros(size(t_hist_l));
    x_vec_hist_h = zeros(size(x_vec_hist_l));
    Phi_mtx_hist_h = zeros(size(Phi_mtx_hist_l));
    i_e_h = zeros(size(i_e_l));
    t_e_h = zeros(size(t_e_l));
    x_e_vec_h = zeros(size(x_e_vec_l));
    Phi_mtx_e_h = zeros(size(Phi_mtx_e_l));
    N = length(t_hist_l);
    Ne = length(t_e_l);

    % Variables that do not change 
    i_e_h = i_e_l;
    t_hist_h = t_hist_l;
    t_e_h = t_e_l;

    % Calculate Transformation matrix T
    T_mtx = [eye(3), zeros(3);
            [0,-1,0; 1, 0, 0; 0, 0, 0], eye(3)];
    T_inv_mtx = [eye(3), zeros(3);
            [0,1,0; -1, 0, 0; 0, 0, 0], eye(3)];


    % Transform coordinates
    for i = 1 : N
        x_vec_h = T_mtx*x_vec_hist_l(i, :)';
        x_vec_hist_h(i, :) = x_vec_h';
    end

    for i = 1 : Ne
        x_e_vec = T_mtx*x_e_vec_l(i, :)';
        x_e_vec_h(i, :) = x_e_vec';
    end

    % Transform the STM
    for i = 1 : N
        Phi_mtx_l = squeeze(Phi_mtx_hist_l(i, :,:));
        Phi_mtx_hist_h(i, :,:) = T_mtx*Phi_mtx_l*T_inv_mtx;
    end

    for i = 1 : Ne
        Phi_mtx_l = squeeze(Phi_mtx_e_l(i, :,:));
        Phi_mtx_e_h(i, :,:) = T_mtx*Phi_mtx_l*T_inv_mtx;
    end

    % Changing the model of the parameters
    params_h = params_l;
    params_h.model.formulation = 'hamiltonian';

end