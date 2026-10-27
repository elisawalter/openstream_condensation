classdef CONDHTM
    %CONDHTM Enumeration of condensation wall heat transfer models
    %
    % This class defines the available :attr:`Inputs.Model.CONDHTM` models
    % for calculating condensation wall heat transfer in the mixture solver
    % when wall heat flux is negative and the local state is two-phase.
    %
    % Models:
    %
    % - NONE      — Disable dedicated condensation HTC model
    % - NUSSELT   — Vertical-film condensation correlation based on film Reynolds number
    % - SHAH      — Shah correlation for condensation in tubes
    % - KIMMUDAWAR — Kim-Mudawar correlation for shear-dominated condensation

    enumeration
        NONE                   % No dedicated condensation HTC model
        NUSSELT                % Nusselt-type vertical film model
        SHAH                   % Shah condensation model
        KIMMUDAWAR             % Kim-Mudawar condensation model
    end
end
