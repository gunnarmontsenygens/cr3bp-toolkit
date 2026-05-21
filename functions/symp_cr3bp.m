function [t_hist, symp_err_hist, symp_err_f, symp_err_max, symp_valid] = symp_cr3bp(t_hist_h, Phi_mtx_hist_h, tol, params_h)
%==========================================================================
%
% Evaluates the symplecticity of the State Transition Matrix (STM) for the
% Circular Restricted Three-Body Problem (CR3BP) in Hamiltonian coordinates.
%
% MODEL DESCRIPTION:
% In Hamiltonian systems, the flow is symplectic, meaning that the STM must
% preserve the canonical symplectic form. For the CR3BP formulated in
% Hamiltonian coordinates, the STM Phi satisfies:
%
%   Phi^T * J * Phi = J
%
% where J is the canonical symplectic matrix.
%
% This property is fundamental to Hamiltonian dynamics and reflects the
% preservation of phase-space structure under time evolution.
%
% This function evaluates the deviation from symplecticity along a trajectory
% by computing the Frobenius norm of:
%
%   E = Phi^T * J * Phi - J
%
% at each time step.
%
% STATE DEFINITION:
% The STM is assumed to correspond to the Hamiltonian state:
%
%   x_vec = [x; y; z; p_x; p_y; p_z]
%
% where (x, y, z) are positions and (p_x, p_y, p_z) are canonical momenta.
%
% SYMPLECTIC FORM:
% The canonical symplectic matrix is defined as:
%
%   J = [ 0   I
%        -I   0 ]
%
% where I is the 3x3 identity matrix.
%
% The symplecticity condition ensures that the STM preserves phase-space
% volume and the geometric structure of the Hamiltonian flow.
%
% OUTPUT DESCRIPTION:
% The function returns:
%   - The time history
%   - The symplectic error at each time step
%   - The final symplectic error
%   - The maximum symplectic error over the trajectory
%   - A logical flag indicating whether the symplectic condition is satisfied
%     within a specified tolerance
%
% NOTES:
% - This check is only valid for Hamiltonian formulations of the CR3BP.
% - The symplectic error should remain close to machine precision for a
%   correctly implemented STM.
% - Large deviations may indicate numerical integration issues or incorrect
%   STM propagation.
%
%
% Author: G. Montseny
% Date: May 5, 2026
%
% INPUT:                 Description                                   Units
%
%  t_hist_h       -   time history (Nx1)                              [-]
%  Phi_mtx_hist_h -   STM history (Nx6x6)                             [-]
%  tol            -   tolerance for symplecticity check               [-]
%  params_h       -   parameter struct (Hamiltonian formulation)      [-]
%
% OUTPUT:                Description                                   Units
%
%  t_hist         -   time history (Nx1)                              [-]
%  symp_err_hist  -   symplectic error history (Nx1)                  [-]
%  symp_err_f     -   final symplectic error                          [-]
%  symp_err_max   -   maximum symplectic error                        [-]
%  symp_valid     -   logical flag (true if error < tol)              [-]
%
%==========================================================================```
    
    % Check if system is given in Hamiltonian coordinates
    switch lower(params_h.model.formulation)
        case 'hamiltonian'
        otherwise
            error('Formulation needs to be Hamiltonian')
    end

    % Preallocation
    t_hist = t_hist_h;
    N = length(t_hist_h);
    symp_err_hist = zeros(N,1);
    
    % Definition of the symplectic identity matrix
    J = [zeros(3), eye(3);
        - eye(3), zeros(3)];

    % Calculate the symplectic error at each time period
    for i = 1 : N
        Phi_mtx = squeeze(Phi_mtx_hist_h(i, :, :));
        symp_err_mtx = Phi_mtx'*J*Phi_mtx - J;
        symp_err = norm(symp_err_mtx,'fro');
        symp_err_hist(i) = symp_err;
    end

    % Calculate remaining outputs
    symp_err_f = symp_err_hist(end);
    symp_err_max = max(symp_err_hist);

    if symp_err_max < tol
        symp_valid = true;
    else
        symp_valid = false;
    end


end