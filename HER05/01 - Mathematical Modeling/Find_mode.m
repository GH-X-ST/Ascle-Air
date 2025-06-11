function modes = Find_mode(A)
% Introduction:
%   Each field contains the natural frequency in rad/s.
%   If a mode cannot be unambiguously identified it is set to NaN.
%
% Author:
%   Hanchen Li (hl3422@ic.ac.uk)
%
% Inputs:
% Outputs:
%   A         - system matrix
%
% Outputs:
%   mode
%             - shortPeriod
%             - dutchRoll
%             - phugoid
%             - rollSubsidence
%             - heave
%             - spiral

%% 1 Eigen–decomposition
lam = eig(A);

% 1.1 Separate real and complex‐conjugate eigenvalues
realRoots     = lam( abs(imag(lam)) < 1e-6 );
complexPairs  = lam( abs(imag(lam)) >= 1e-6 );

% 1.2 Keep only one representative of each conjugate pair
complexPairs  = complexPairs( imag(complexPairs) > 0 );

% 1.3 Natural frequency and damping for complex pairs
wnComplex = abs(complexPairs);      % ω_n = |λ|

% 1.4 Sort descending by ω_n
[wnSorted,idx] = sort(wnComplex,'descend');
sortedPairs    = complexPairs(idx);

%% 2 Complex modes
modes.shortPeriod = pick(wnSorted,1);
modes.dutchRoll   = pick(wnSorted,2);
modes.phugoid     = pick(wnSorted,numel(wnSorted));

%% 3 Real modes
% 3.1 sort real roots by their (negative) real part (most negative last)
[~,iReal]  = sort(real(realRoots));   % ascending is increasingly negative
sortedReal = realRoots(iReal);

% 3.2 roll-subsidence is fastest (most negative) real root
modes.rollSubsidence = abs(sortedReal(1));

% 3.3 heave (w-subsidence) is second-most-negative real root
modes.heave = pick(abs(sortedReal),2);

% 3.4 spiral is real root closest to the imaginary axis
nzReal       = realRoots( abs(real(realRoots)) > 1e-4 );   % tolerance
[~,iS]       = min(abs(nzReal));
modes.spiral = abs(real(nzReal(iS)));

%% A Helper
% safe element picker (returns NaN if index is out of bounds)
    function val = pick(arr,k)
        if k >= 1 && k <= numel(arr)
            val = arr(k);
        else
            val = NaN;
        end
    end

end