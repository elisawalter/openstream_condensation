classdef CondensationHTCUnitTest < matlab.unittest.TestCase
    % CONDENSATIONHTCUNITTEST
    %
    % Unit tests for the condensation heat transfer correlations in
    % Solvers.Mixture.Mixture.CONDHTC.
    %
    % The class verifies that:
    %
    %   1. All CONDHTM models give finite, positive coefficients over the
    %      two-phase quality range.
    %   2. SHAH reproduces the Shah (1979) correlation.
    %   3. NUSSELTSHAH is the larger of NUSSELT and SHAH.
    %   4. KIMMUDAWAR includes the two-phase multiplier, i.e. it exceeds
    %      the bare film term 0.048 Re_f^0.69 Pr_l^0.34 k_l/D.
    %   5. The coefficients are expanded to all walls.
    %   6. KIMMUDAWAR is continuous across the laminar/turbulent and
    %      annular/slug transitions, and finite as x -> 0.
    %
    % This class does not construct or run any OpenSTREAM solver.

    properties (SetAccess = private)

        x         % Vapor quality [-]
        Ref       % Film Reynolds number [-]
        G         % Mass flux [kg/m^2/s]
        D         % Hydraulic diameter [m]
        g         % Gravitational acceleration [m/s^2]
        props     % Saturation properties of water at 4 bar

    end

    methods (TestClassSetup)

        function configureTestEnvironment(testCase)
            % CONFIGURETESTENVIRONMENT
            %
            % Add the OpenSTREAM repository to the MATLAB path and define
            % PICON-V-like conditions: water at 4 bar, 40 mm tube, G = 30 kg/m^2/s.

            openstreamRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(openstreamRoot));

            testCase.x   = (0.05:0.1:0.95).';
            testCase.G   = 30*ones(size(testCase.x));
            testCase.D   = 0.04;
            testCase.g   = 9.81;
            testCase.props = struct( ...
                'RHOL',  922.0, ...
                'RHOV',  2.16, ...
                'MUL',   1.85E-4, ...
                'MUV',   1.37E-5, ...
                'KL',    0.686, ...
                'PRL',   1.16, ...
                'SIGMA', 0.0489, ...
                'PR',    4E5/22.064E6);
            testCase.Ref = testCase.G.*(1-testCase.x).*testCase.D./testCase.props.MUL;
        end

    end

    methods (Test)

        function finitePositive(testCase)
            % FINITEPOSITIVE
            %
            % Verify finite, positive coefficients for all models except NONE.

            for m = setdiff(string(enumeration('InputEnums.CONDHTM')).', "NONE")
                h = testCase.evaluate(m);
                testCase.verifyTrue(all(isfinite(h) & h > 0), ...
                    sprintf("CONDHTM %s gives non-finite or non-positive HTC.", m));
            end

            testCase.verifyEqual(testCase.evaluate("NONE"), zeros(size(testCase.x)), ...
                "CONDHTM NONE must give zero HTC.");
        end

        function shah1979(testCase)
            % SHAH1979
            %
            % Verify SHAH against the Shah (1979) correlation written out.

            p = testCase.props;
            hlo = 0.023.*(testCase.G.*testCase.D./p.MUL).^0.8.*p.PRL.^0.4.*p.KL./testCase.D;
            hRef = hlo.*((1-testCase.x).^0.8 + 3.8.*testCase.x.^0.76.*(1-testCase.x).^0.04./p.PR.^0.38);

            testCase.verifyEqual(testCase.evaluate("SHAH"), hRef, 'RelTol', 1E-12, ...
                "SHAH does not reproduce Shah (1979).");
        end

        function nusseltShahIsMaximum(testCase)
            % NUSSELTSHAHISMAXIMUM
            %
            % Verify NUSSELTSHAH = max(NUSSELT, SHAH).

            testCase.verifyEqual(testCase.evaluate("NUSSELTSHAH"), ...
                max(testCase.evaluate("NUSSELT"), testCase.evaluate("SHAH")), ...
                "NUSSELTSHAH must be the larger of NUSSELT and SHAH.");
        end

        function kimMudawarMultiplier(testCase)
            % KIMMUDAWARMULTIPLIER
            %
            % Verify that KIMMUDAWAR includes the phi_g/X_tt enhancement.

            p = testCase.props;
            hFilm = 0.048.*testCase.Ref.^0.69.*p.PRL.^0.34.*p.KL./testCase.D;

            testCase.verifyGreaterThan(testCase.evaluate("KIMMUDAWAR")./hFilm, 1, ...
                "KIMMUDAWAR is missing the two-phase multiplier.");
        end

        function kimMudawarContinuous(testCase)
            % KIMMUDAWARCONTINUOUS
            %
            % Verify that KIMMUDAWAR has no jump at the laminar/turbulent
            % transitions and the annular/slug regime change. With
            % G = 10 kg/m^2/s, Re_g passes 2000-4000 at x = 0.06-0.12; a
            % step of 1E-4 in x must change h by < 1% over the quality range.

            G  = 10;
            xf = (0.01:1E-4:0.99).';
            Rf = G.*(1-xf).*testCase.D./testCase.props.MUL;
            h  = Solvers.Mixture.Mixture.CONDHTC("KIMMUDAWAR", xf, Rf, G, ...
                testCase.D, testCase.g, testCase.props);

            testCase.verifyLessThan(max(abs(diff(h))./h(1:end-1)), 0.01, ...
                "KIMMUDAWAR jumps at a regime transition.");

            % Slug/bubbly term keeps the HTC finite at vanishing quality.
            h0 = Solvers.Mixture.Mixture.CONDHTC("KIMMUDAWAR", 1E-6, ...
                G*testCase.D/testCase.props.MUL, G, testCase.D, testCase.g, testCase.props);
            testCase.verifyGreaterThan(h0, 500, ...
                "KIMMUDAWAR vanishes at low quality (slug/bubbly term missing).");
        end

        function expandToWalls(testCase)
            % EXPANDTOWALLS
            %
            % Verify that the coefficient takes the size of the per-wall
            % film Reynolds number.

            Ref2 = [testCase.Ref, 0.5*testCase.Ref];
            for m = string(enumeration('InputEnums.CONDHTM')).'
                h = Solvers.Mixture.Mixture.CONDHTC(m, testCase.x, Ref2, ...
                    testCase.G, testCase.D, testCase.g, testCase.props);
                testCase.verifySize(h, size(Ref2), ...
                    sprintf("CONDHTM %s does not expand to all walls.", m));
            end
        end

    end

    methods (Access = private)

        function h = evaluate(testCase, condhtm)
            % EVALUATE Condensation HTC for the test conditions [W/m^2/K]

            h = Solvers.Mixture.Mixture.CONDHTC(condhtm, testCase.x, testCase.Ref, ...
                testCase.G, testCase.D, testCase.g, testCase.props);
        end

    end

end
