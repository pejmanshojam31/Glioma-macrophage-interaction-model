function [prcc] = TimeDependentPRCC(M, N, A, OutputOfInterest, timePoints, fixedSpatialIdx)
    % Initialize matrices
    K = zeros(N, M + 1);  % Include output column
    R = zeros(N, M + 1);  % For ranks
    prcc = zeros(M, length(timePoints));  % Store PRCC values

    % Loop over time points
    for t = 1:length(timePoints)
        % Assign parameter values and corresponding output at fixed spatial location
        K(:, 1:M) = A(:, 1:M);
        for n = 1:N
            K(n, M + 1) = OutputOfInterest(n, fixedSpatialIdx, timePoints(t));
        end

        % Rank the data
        for m = 1:M + 1
            [~, i] = sort(K(:, m));
            r = 1:N;
            R(:, m) = r(i);
        end

        % Compute correlation matrix and PRCC
        C = corrcoef(R);
        if det(C) <= 10^-16
            B = pinv(C);
            fprintf('C is nearly singular at time point %d.\n', t);
        else
            B = inv(C);
        end

        % Calculate PRCCs
        for w = 1:M
            prcc(w, t) = (-B(w, M + 1)) / sqrt(B(w, w) * B(M + 1, M + 1));
        end
    end
end