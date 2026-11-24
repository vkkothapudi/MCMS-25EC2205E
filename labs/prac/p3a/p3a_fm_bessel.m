%% P3A: FM modulation - Bessel spectrum, Carson's rule and the discriminator (MCMS 25EC2205E)
% Generates tone FM at several modulation indices, measures the spectral line heights,
% the 98 % power bandwidth and the deviation recovered by a quadrature discriminator.
% Core MATLAB only (no toolboxes needed).
clear; clc;
%% 1. Parameters: a 1 kHz tone, complex baseband FM (carrier at 0 Hz)
fs = 200e3;                    % sample rate (S/s)
fm = 1e3;                      % message tone (Hz)
N = fs;                        % 1 s of samples: 1 Hz FFT bins, every line falls on a bin
t = (0:N-1)'/fs;
betas = [0.5 1.0 2.405 5.0];
%% 2. Measure every beta
fprintf('%6s %8s %8s %8s %8s %8s %9s %7s %8s\n', 'beta', 'df(kHz)', 'J0 meas', 'J0 th', ...
        'J1 meas', 'J1 th', 'B98(kHz)', 'Carson', 'df meas');
for b = betas
    s = fm_signal(b, fm, t);
    a = lines(s, fm, N, 12); c = 13;           % index of the carrier line
    fi = discriminator(s, fs);
    fprintf('%6.3f %8.3f %8.4f %8.4f %8.4f %8.4f %9.1f %7.1f %8.3f\n', b, b*fm/1e3, a(c), ...
            abs(besselj(0, b)), a(c+1), abs(besselj(1, b)), bw98(a, fm)/1e3, 2*(b*fm + fm)/1e3, max(fi)/1e3);
end
%% 3. Carrier-level sweep: the carrier line follows |J0(beta)|
bsweep = 0:0.25:6; j0m = zeros(size(bsweep));
for k = 1:numel(bsweep)
    j0m(k) = lines(fm_signal(bsweep(k), fm, t), fm, N, 0);
end
fprintf('Carrier sweep beta 0..6: max |J0 measured - J0 theory| = %.1e\n', max(abs(j0m - abs(besselj(0, bsweep)))));
fprintf('First carrier null at beta = 2.405: carrier = %.4f\n', lines(fm_signal(2.405, fm, t), fm, N, 0));
%% 4. Figure
figure('Position', [100 100 1100 650]);
cols = {'b', 'm'}; bs = [1.0 5.0]; n = -12:12;
for k = 1:2
    subplot(2,2,k); a = lines(fm_signal(bs(k), fm, t), fm, N, 12);
    stem(n*fm/1e3, a, cols{k}, 'filled'); grid on; hold on;
    B = (bs(k)*fm + fm)/1e3; plot([-B -B B B], [0 1 1 0], '--', 'Color', [1 0.5 0]);
    title(sprintf('Spectrum, beta = %.1f (dashed: Carson band)', bs(k)));
    xlabel('Frequency offset (kHz)'); ylabel('Line amplitude');
end
subplot(2,2,3); bb = linspace(0, 6, 400);
plot(bb, abs(besselj(0, bb)), 'k-', bsweep, j0m, 'ro'); grid on;
xlabel('beta'); ylabel('Carrier amplitude'); title('Carrier line vs beta'); legend('|J0(beta)| theory', 'measured carrier');
subplot(2,2,4); fi = discriminator(fm_signal(5.0, fm, t), fs);
plot(t(2:2001)*1e3, fi(1:2000)/1e3, 'Color', [0 0.6 0.2]); grid on;
xlabel('Time (ms)'); ylabel('Instantaneous frequency (kHz)'); title('Discriminator output, beta = 5 (df = 5 kHz)');

%% Local functions (must stay at the end of the script, MATLAB R2016b+)
function s = fm_signal(beta, fm, t)
    s = exp(1j*beta*sin(2*pi*fm*t));        % instantaneous frequency = beta*fm*cos(2*pi*fm*t)
end
function a = lines(s, fm, N, nmax)
    S = fft(s)/N;                           % amplitude of the line at n*fm, n = -nmax..nmax
    k = mod((-nmax:nmax)*fm, N) + 1;
    a = abs(S(k)).';
end
function B = bw98(a, fm)
    p = a.^2; c = (numel(a)+1)/2; B = NaN;  % smallest +-n*fm band holding 98 % of the power
    for n = 0:c-1
        if sum(p(c-n:c+n)) >= 0.98*sum(p), B = 2*n*fm; return; end
    end
end
function fi = discriminator(s, fs)
    fi = angle(s(2:end).*conj(s(1:end-1)))*fs/(2*pi);   % instantaneous frequency (Hz)
end
