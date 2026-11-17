classdef DriftFluxVoidUnitTest < matlab.unittest.TestCase
    % DRIFTFLUXVOIDUNITTEST
    %
    % Unit tests for the bounded drift-flux void fraction
    % Solvers.Mixture.Mixture.VFDRIFTBOUNDED used by the BESTION void model.
    %
    % The class verifies that:
    %
    %   1. The void fraction stays in [0,1] for upflow and downflow drift
    %      velocities, also for annular flow with a thin liquid film.
    %   2. The drift-flux relation vf (C0 j + ugj0 (1-vf)) = JG is solved.
    %   3. With C0 = 1 the phase slip velocity ug - ul equals ugj0.
    %   4. Zero vapor flux gives zero void fraction.
    %
    % This class does not construct or run any OpenSTREAM solver.

    methods (TestClassSetup)

        function configureTestEnvironment(testCase)
            % CONFIGURETESTENVIRONMENT
            %
            % Add the OpenSTREAM repository to the MATLAB path.

            openstreamRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(openstreamRoot));
        end

    end

    methods (Test)

        function boundedUpAndDownflow(testCase)
            % BOUNDEDUPANDDOWNFLOW
            %
            % Verify 0 <= vf <= 1 for annular condensation conditions
            % (jg ~ 20 m/s, thin film) with Bestion-like drift velocities of
            % both signs, and for bubbly conditions.

            JG = [24 22 21 0.1 0.5].';
            JL = [1.5E-4 3E-3 5E-3 1 0.5].';
            for ugj0 = [-2.8 2.8]
                vf = Solvers.Mixture.Mixture.VFDRIFTBOUNDED(JG, JL, 1, ugj0);
                testCase.verifyGreaterThanOrEqual(vf, 0);
                testCase.verifyLessThanOrEqual(vf, 1, ...
                    sprintf("Void fraction exceeds 1 for ugj0 = %g m/s.", ugj0));
            end
        end

        function solvesDriftFluxRelation(testCase)
            % SOLVESDRIFTFLUXRELATION
            %
            % Verify the residual of the drift-flux relation.

            JG = [24 0.1 0.5].'; JL = [3E-3 1 0.5].'; C0 = 1.1; ugj0 = [-2.8; 0.3; 0.3];
            vf = Solvers.Mixture.Mixture.VFDRIFTBOUNDED(JG, JL, C0, ugj0);
            res = vf.*(C0.*(JG+JL) + ugj0.*(1-vf)) - JG;

            testCase.verifyLessThan(abs(res), 1E-9);
        end

        function boundedSlip(testCase)
            % BOUNDEDSLIP
            %
            % Verify ug - ul = ugj0 for C0 = 1, i.e. the liquid velocity stays
            % finite as vf -> 1 in downflow.

            JG = [24 21].'; JL = [1.5E-4 5E-3].'; ugj0 = -2.8;
            vf = Solvers.Mixture.Mixture.VFDRIFTBOUNDED(JG, JL, 1, ugj0);
            slip = JG./vf - JL./(1-vf);

            testCase.verifyEqual(slip, ugj0*ones(size(JG)), 'AbsTol', 1E-3);
        end

        function noVapor(testCase)
            % NOVAPOR
            %
            % Verify vf = 0 without vapor flux.

            vf = Solvers.Mixture.Mixture.VFDRIFTBOUNDED(0, 1, 1, 0.3);
            testCase.verifyEqual(vf, 0);
        end

    end

end
