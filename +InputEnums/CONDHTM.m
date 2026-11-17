classdef CONDHTM
    %CONDHTM Enumeration of condensation wall heat transfer models
    %
    % This class defines the available :attr:`Inputs.Model.CONDHTM` models
    % for calculating condensation wall heat transfer in the mixture solver
    % when wall heat flux is negative and the local state is two-phase.
    % The correlations are evaluated by
    % :meth:`Solvers.Mixture.Mixture.CONDHTC`.
    %
    % Models:
    %
    % - NONE        — Disable dedicated condensation HTC model
    % - NUSSELT     — Gravity-driven vertical-film condensation (Nusselt laminar, Kutateladze wavy, Labuntsov turbulent) based on film Reynolds number
    % - SHAH        — Shah correlation for film condensation inside tubes (:cite:t:`shah1979`)
    % - SHAH2009    — Shah correlation for vertical tubes with gravity- and shear-dominated regimes (:cite:t:`shah2009`)
    % - KIMMUDAWAR  — Kim-Mudawar annular-flow correlation for shear-dominated condensation (:cite:t:`kim2013`), with the two-phase multiplier of :cite:t:`kim2012`
    % - NUSSELTSHAH — Larger of NUSSELT and SHAH, which adds vapor-shear enhancement at high mass flux to the gravity-driven film model

    enumeration
        NONE                   % No dedicated condensation HTC model
        NUSSELT                % Nusselt-type vertical film model
        SHAH                   % Shah (1979) condensation model
        SHAH2009               % Shah (2009) vertical tube condensation model
        KIMMUDAWAR             % Kim-Mudawar (2013) annular condensation model
        NUSSELTSHAH            % Larger of Nusselt and Shah (1979)
    end
end
