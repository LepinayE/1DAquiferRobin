function [U] = flowVelocity(t_vec, period)

% flowVelocity - Creates a vector of the non-dimentional flow velocity for
% the injection and extraction fluid.
%
% Syntax:  [U] = flowVelocity(t_vec,P)
%
% Inputs:
%    t_vec - Discretised non-dim time vector over which full system is run
%    P - Period of time for 1 oscillation (years/ number of total cycles)
% 
% Outputs:
%    U - Flow velocity in aquifer at each time step k
%
% Runs in : stepaquifer.m, solveAquifertemp.m

%--------------------------------------------------------------------------
% Author: Emma Lepinay
% Email: el547@cam.ac.uk
% Date: 31/08/2022; Last revision: 
% Version: R2022a

%------------- BEGIN CODE -------------------------------------------------
    U = zeros(size(t_vec));
    t = t_vec;
    t2 = t-period.*floor(t ./ period); % time during each period

    % Velocity is a step function for injection and extraction
    U = 1 .* (t2< period/2)  - 1 .* (t2>=period/2);

end
%------------- END OF CODE ------------------------------------------------

