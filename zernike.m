% (C) 2025, Muthu Annamalai, <muthu@panmo.cloud>
% This file is released under Apache 2.0 license
% If you find this program useful in your research work please cite the article:
% M. Annamalai, "Techniques to Improve Computation of Zernike Polynomials,"
% Frontiers in Optics - Laser Science, Denver, CO (2025).

%sample script
1+1;

% ---- Zernike utilities ----

function R = radial(n, m, rho)
% Compute the Zernike radial polynomial R_n^m(rho) for each value in rho.
% Parameters:
%   n >= 0, |m| <= n, (n-m) even
% rho: array (same size for output)
    if mod(n - m, 2) ~= 0
        R = zeros(size(rho));
        return;
    end
    R = zeros(size(rho));
    k_max = (n - m) / 2;
    for k = 0:k_max
        num = (-1)^k * factorial(n - k);
        den = factorial(k);
        invTerm = 1.0 / ( factorial((n + m)/2 - k) * factorial((n - m)/2 - k) );
        R = R + (num/den) * invTerm .* rho.^(n - 2*k);
    end
end

function field = zernike(n, m, rho, theta, fieldcache)
% Compute the full Zernike polynomial Z_n^m(rho, theta).
% Optional simple cache & 3-term recurrence when available.
    if nargin < 5
        fieldcache = [];
    end

    if ~isempty(fieldcache)
        key = cacheKey(n, m);
        if isKey(fieldcache, key)
            field = fieldcache(key);
            return;
        end

        k1 = cacheKey(n-1, m-1);
        k2 = cacheKey(n-1, m+1);
        k3 = cacheKey(n-2, m);
        if isKey(fieldcache, k1) && isKey(fieldcache, k2) && isKey(fieldcache, k3)
            t1 = fieldcache(k1);
            t2 = fieldcache(k2);
            t3 = fieldcache(k3);
            field = rho .* (t1 + t2) - t3;
            fieldcache(key) = field;
            return;
        end
    end

    R = radial(n, abs(m), rho);
    field = R .* sin( (m >= 0) .* (pi/2 - abs(m).*theta) + (m < 0) .* (abs(m).*theta) );

    if ~isempty(fieldcache)
        fieldcache(cacheKey(n, m)) = field;
    end
end

function key = cacheKey(n, m)
    key = sprintf('%d_%d', n, m);
end

function basis = generate_basis(n_max, resolution)
% Generate all Zernike basis modes up to order n_max on a square grid.
% Returns a cell array of size K with each entry resolution x resolution.
    persistent global_fieldcache
    if isempty(global_fieldcache)
        global_fieldcache = containers.Map('KeyType','char','ValueType','any');
    end

    % normalized grid from -1 to 1 (Y flipped so positive up)
    x = linspace(-1, 1, resolution);
    y = linspace(-1, 1, resolution);
    [X, Y] = meshgrid(x, -y);
    rho = sqrt(X.^2 + Y.^2);
    theta = atan2(Y, X);
    mask = (rho <= 1.0);

    % Optional recursion shortcut like Python version
    if n_max > 10
        basis = generate_basis(max(0, n_max - 2), resolution);
        range_n = (n_max-2):n_max;
    else
        basis = {};
        range_n = 0:n_max;
    end

    for n = range_n
        for m = -n:n
            if mod(n - m, 2) == 0
                Z = zernike(n, m, rho, theta, global_fieldcache);
                Z(~mask) = 0.0;
                basis{end+1} = Z; %#ok<AGROW>
            end
        end
    end
end

function coeffs = compute_moments(image, basis)
% Least-squares fit of image onto provided Zernike basis.
    N = numel(image);
    K = numel(basis);
    A = zeros(N, K);
    b = image(:);
    for k = 1:K
        A(:, k) = basis{k}(:);
    end
    % Equivalent to numpy.linalg.lstsq: use backslash; falls back to LS
    coeffs = A \ b;
end

function img = recreate_image(coeffs, basis)
% Reconstruct image from Zernike coefficients and basis.
    img = zeros(size(basis{1}));
    for k = 1:numel(coeffs)
        img = img + coeffs(k) .* basis{k};
    end
end

function M = circular_mask(resolution)
% Return a circular mask (1 inside unit circle, 0 outside).
    x = linspace(-1, 1, resolution);
    y = linspace(-1, 1, resolution);
    [X, Y] = meshgrid(x, -y);
    rho = sqrt(X.^2 + Y.^2);
    M = double(rho <= 1.0);
end

% Zernike Demo:
    % Example usage (equivalent to Python __main__)
    N = 256;                               % grid size
    for n_max = 10:2:16                    % 10, 12, 14, 16
        t0 = tic;
        basis = generate_basis(n_max, N);
        fprintf('Elapsed for %d = %.3f s\n', n_max, toc(t0));

        % 2) Test image: unit-circle mask
        img = circular_mask(N);

        % 3) Compute Zernike moments (least squares)
        coeffs = compute_moments(img, basis);
        fprintf('Computed %d coefficients.\n', numel(coeffs));

        % 4) Reconstruct image
        img_rec = recreate_image(coeffs, basis);

        % 5) Relative reconstruction error
        err = norm(img(:) - img_rec(:)) / norm(img(:));
        fprintf('@Order %d Relative reconstruction error: %.3e\n', n_max, err);
    end
