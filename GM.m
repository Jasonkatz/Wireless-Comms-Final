function [ Y ] = GM( x, mu, omega )

    Y = 1 / det(pi * omega) * exp(-(x - mu)' * pinv(omega) * (x - mu));

end

