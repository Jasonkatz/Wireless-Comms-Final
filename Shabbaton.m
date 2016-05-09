% ECE-408 Final Project
% Elie Lerea and Jason Katz - Team Shabbaton
% Massive MIMO Channel Estimator Simulation

%% Simulation Setup
clc, clear all, clear global, close all;
dbstop if error;

msgM = 2; % BPSK
k = log2(msgM);
B = 4; % Cells
K = 80; % User environments per cell
N = 16; % Antennas per BS
T = 80; % Number of pilot symbols
xi = 1 / T;

%% Data Generation

% Generate S matrix (TxK, line up B times -> TxBK)
S = [];
for i = 1:(B*K)
    bits = randi([0 1], k * T, 1); % Generate random bits, pass these out of function, unchanged
    syms = bi2de(reshape(bits,k,length(bits)/k).','left-msb')';
    msg = qammod(syms, msgM)';
    S = [S msg];
end
S = S * sqrt(xi);

% Generate H matrix (KxN, stack up B times -> BKxN)
H = [];
for i = 1:B
    H = [H ; randn(K, N)];
end

% Take DFT of channels
F = dftmtx(N);
H_ = H * F;

Y = S * H;
Y_ = Y * F;

% This is all for testing
a = rand(320, 256);
v = rand(320, 256);
V = ones(K, 1); % V0 is 1; THIS MIGHT BE T NOT K CAREFUL
w = zeros(K, N); % w0 is 0; THIS MIGHT BE T NOT K CAREFUL
sigSq = zeros(K, N);
R = zeros(K, N);
Q = zeros(B*K, N);
epsilon = 1e-6;

L = 4;
etaRo = .5 * ones(1, L);
etaSigma = 2 .^ [0:(L-1)];
etaRo = double(etaRo);
etaSigma = double(etaSigma);

for n = 1:N
    
    for t = 1:20 % 20 iterations

        Vnew = repmat(xi * sum(v(:, n)), K, 1);
        wnew = S * a(:, n) - Vnew ./ V .* (Y_(:, n) - w(:, n)); % NOT SURE ABOUT THIS BUT WE DONT CARE
        sigSq = repmat(Vnew / (xi * T), B, 1);

        temp = zeros(B*K, 1);
        for k = 1:B*K
            for i = 1:T % THIS MIGHT BE K NOT T CAREFUL
                temp(k) = temp(k) + (S(i, k) * (Y_(i, n) - wnew(i)));
            end
        end
        R = a(:, n) + temp / (xi * T);
        
        anew = [];
        vnew = [];
        for k = 1:(B*K)
            numA = 0;
            denomA = 0;
            numB = 0;
            denomB = 0;
            for l = 1:L
                % WHAT THE FUCK
%                 numA = numA + etaRo(l) * GM(0, R(k), sigSq(k) + etaSigma(l)) * R(k) * etaSigma(l) / (sigSq(k) + etaSigma(l));
%                 denomA = denomA + etaRo(l) * GM(0, R(k), sigSq(k) + etaSigma(l));
%                 numB = numB + etaRo(l) * GM(0, R(k), sigSq(k) + etaSigma(l)) * (etaSigma(l) * sigSq(k) * (sigSq(k) + etaSigma(l)) + abs(R(k))^2*etaSigma(l)^2) / (sigSq(k) + etaSigma(l));
%                 denomB = denomB + etaRo(l) * GM(0, R(k), sigSq(k) + etaSigma(l));
                numA = R(k) * etaSigma(l) / (sigSq(k) + etaSigma(l));
                denomA = 1;
                numB = (etaSigma(l) * sigSq(k) * (sigSq(k) + etaSigma(l)) + abs(R(k))^2*etaSigma(l)^2) / (sigSq(k) + etaSigma(l));
                denomB = 1;
            end
            anew = [anew ; numA / denomA];
            vnew = [vnew ; numB / denomB];
        end
        anew(isnan(anew)) = 0;
        vnew(isnan(vnew)) = 0;
        
        aold = a;
        V = Vnew;
        w(:, n) = wnew;
        a(:, n) = anew;
        v(:, n) = vnew;
        
        for k = 1:(B*K)
            Q(k, n) = GM(H_(k, n), a(k), v(k));
        end
        
        for l = 1:L
            sumFk = 0;
            sumWhateva = 0;
            for k = 1:(B*K)
                numFk = etaRo(l) * GM(0, a(k), etaSigma(l) + v(k));
                denomFk = 0;
                gamma = etaSigma(l) / (etaSigma(l) + v(k)) * a(k);
                zeta = etaSigma(l) * v(k) / (etaSigma(l) + v(k));
                for lp = 1:L
                    denomFk = denomFk + etaRo(lp) * GM(0, a(k), etaSigma(lp) + v(k));
                end
                sumFk = sumFk + numFk / denomFk;
                sumWhateva = sumWhateva + numFk / denomFk * (abs(gamma)^2 + zeta);
            end
            etaRo(l) = 1 / (B*K) * sumFk;
            etaSigma(l) = sumWhateva / (B * T * etaRo(l));
        end

        % Compute error
        err = double(immse(aold(:, n), a(:, n))) * 320
        if err < epsilon
            break;
        end
        
    end

end