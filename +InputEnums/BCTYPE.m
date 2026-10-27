classdef BCTYPE
    %BCTYPE Enumeration of wall thermal boundary-condition types
    %
    % This class defines available wall thermal boundary-condition types
    % for the mixture solver.
    %
    % Types:
    %
    % - HFLUX    — Prescribed wall heat flux profile
    % - TWALLIN  — Prescribed inner-wall temperature profile
    % - TWALLOUT — Prescribed outer-wall temperature profile

    enumeration
        HFLUX
        TWALLIN
        TWALLOUT
    end
end
