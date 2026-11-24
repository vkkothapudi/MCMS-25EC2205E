%% P3B: Phase modulation - PM vs FM, the narrowband phase modulator and PM demodulation (MCMS 25EC2205E)
% Compares the frequency deviation of PM and FM as the tone frequency changes, shows what a
% triangle message does to each, measures the error of the narrowband (Armstrong) phase modulator
% and recovers the message from the phase. Core MATLAB only (no toolboxes needed).
clear; clc;
%% 1. Parameters
fs = 200e3;                      % sample rate (S/s)
N = fs;                          % 1 s of samples
t = (0:N-1)'/fs;
kp = 2.0;                        % PM phase sensitivity (rad per volt)
kf = 5e3;                        % FM frequency sensitivity (Hz per volt)
Am = 1.0;                        % message amplitude (V)
tone = @(f) Am*cos(2*pi*f*t);
tri  = @(f) Am*(2*abs(2*mod(f*t + 0.25, 1) - 1) - 1);     % triangle, +-Am
pm   = @(m) exp(1j*kp*m);                                % phase follows m(t)
fmod = @(m) exp(1j*2*pi*kf*cumsum(m)/fs);                % frequency follows m(t)
%% 2. Peak deviation vs tone frequency: FM fixed, PM grows with fm
fprintf('Part 1: tone message, Am = 1 V, kp = 2 rad/V, kf = 5 kHz/V\n');
fprintf('%7s %10s %8s %10s %8s %7s %7s\n', 'fm(Hz)', 'dfFM meas', 'dfFM th', 'dfPM meas', 'dfPM th', 'betaFM', 'betaPM');
fms = [250 500 1000 2000 4000]; dfm = zeros(size(fms)); dpm = dfm;
for k = 1:numel(fms)
    f = fms(k); m = tone(f);
    dfm(k) = max(inst_freq(fmod(m), fs)); dpm(k) = max(inst_freq(pm(m), fs));
    fprintf('%7d %10.1f %8.1f %10.1f %8.1f %7.2f %7.2f\n', f, dfm(k), kf*Am, dpm(k), kp*Am*f, kf*Am/f, kp*Am);
end
%% 3. Triangle message: PM frequency jumps between two levels, FM sweeps smoothly
ftri = 500;
fp = inst_freq(pm(tri(ftri)), fs); ff = inst_freq(fmod(tri(ftri)), fs);
th = kp*4*Am*ftri/(2*pi);                     % |dm/dt| = 4*Am*f for a triangle
mid = abs(fp) > 0.5*th;                       % samples away from the jumps
fprintf('\nPart 2: triangle message, %d Hz\n', ftri);
fprintf('PM: frequency jumps between +-%.1f Hz (theory +-%.1f Hz), %.1f %% of samples within 1 %% of a level\n', ...
        median(abs(fp(mid))), th, mean(abs(abs(fp) - th) < 0.01*th)*100);
fprintf('FM: frequency sweeps smoothly between %.0f and %.0f Hz (theory +-%.0f Hz)\n', min(ff), max(ff), kf*Am);
%% 4. Narrowband phase modulator: cos(wc t) - theta*sin(wc t) at baseband = 1 + j*theta
fprintf('\nPart 3: narrowband phase modulator, s = 1 + j*theta (theta = k*m, tone 1 kHz)\n');
fprintf('%8s %9s %8s %9s %8s %7s\n', 'theta pk', 'phase pk', 'error %', 'th err %', 'env max', 'th env');
nb_k = [0.1 0.2 0.5 1.0]; nb_err = zeros(size(nb_k));
for k = 1:numel(nb_k)
    s = 1 + 1j*nb_k(k)*tone(1e3);
    ph = max(angle(s)); nb_err(k) = (nb_k(k) - ph)/nb_k(k)*100;
    fprintf('%8.2f %9.4f %8.2f %9.2f %8.4f %7.4f\n', nb_k(k), ph, nb_err(k), ...
            (nb_k(k) - atan(nb_k(k)))/nb_k(k)*100, max(abs(s)), sqrt(1 + nb_k(k)^2));
end
%% 5. PM demodulator: the phase of s(t) divided by kp gives m(t) back
sine = Am*sin(2*pi*1e3*t);                    % a sine starts at 0, so the phase starts unwrapped
m = sine*1.4;  mh = angle(pm(m))/kp;          % peak phase 2.8 rad, still inside +-pi
fprintf('\nPart 4: PM demodulation, peak phase %.1f rad: rms error = %.1e V\n', kp*1.4, sqrt(mean((mh - m).^2)));
m2 = sine*2.0; mw = angle(pm(m2))/kp; mh2 = unwrap(angle(pm(m2)))/kp;   % peak phase 4 rad
fprintf('Peak phase %.1f rad (beyond pi): wrapped error = %.2f V, after unwrap = %.1e V\n', kp*2.0, ...
        sqrt(mean((mw - m2).^2)), sqrt(mean((mh2 - m2).^2)));
%% 6. Figure
figure('Position', [100 100 1100 650]);
subplot(2,2,1); tm = t(2:801)*1e3;
plot(tm, ff(1:800)/1e3, 'r', tm, fp(1:800)/1e3, 'm'); grid on;
xlabel('Time (ms)'); ylabel('Instantaneous frequency (kHz)'); title('Triangle message, 500 Hz'); legend('FM of triangle', 'PM of triangle');
subplot(2,2,2);
plot(fms, dfm/1e3, 'ro-', fms, dpm/1e3, 'ms-', fms, kp*Am*fms/1e3, 'k--'); grid on;
xlabel('Tone frequency fm (Hz)'); ylabel('Peak deviation (kHz)'); title('Deviation vs tone frequency');
legend('FM (measured)', 'PM (measured)', 'PM theory kp*Am*fm', 'Location', 'northwest');
subplot(2,2,3); kk = linspace(0.02, 1.2, 100);
plot(kk, (kk - atan(kk))./kk*100, 'k-', nb_k, nb_err, 'go'); grid on;
xlabel('Peak phase theta (rad)'); ylabel('Phase error (%)'); title('Narrowband phase modulator error'); legend('theory', 'measured', 'Location', 'northwest');
subplot(2,2,4); tt = t(1:600)*1e3;
plot(tt, m2(1:600), 'k--', tt, mw(1:600), 'Color', [1 0.5 0]); hold on; plot(tt, mh2(1:600), 'b'); grid on;
xlabel('Time (ms)'); ylabel('Message (V)'); title('PM demodulation, peak phase 4 rad');
legend('message sent', 'phase / kp (wrapped)', 'unwrapped phase / kp', 'Location', 'southwest');

%% Local functions (must stay at the end of the script, MATLAB R2016b+)
function fi = inst_freq(s, fs)
    fi = angle(s(2:end).*conj(s(1:end-1)))*fs/(2*pi);   % instantaneous frequency (Hz)
end
