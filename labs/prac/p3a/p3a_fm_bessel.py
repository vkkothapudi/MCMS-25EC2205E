"""P3A: FM modulation - Bessel spectrum, Carson's rule and the discriminator (MCMS 25EC2205E) - Python version.
Generates tone FM at several modulation indices, measures the spectral line heights,
the 98 % power bandwidth and the deviation recovered by a quadrature discriminator.
Needs: numpy, scipy, matplotlib."""
import numpy as np
from scipy.special import jv
import matplotlib.pyplot as plt

# 1. Parameters: a 1 kHz tone, complex baseband FM (carrier at 0 Hz)
fs = 200e3                     # sample rate (S/s)
fm = 1e3                       # message tone (Hz)
N = int(fs)                    # 1 s of samples: 1 Hz FFT bins, every line falls on a bin
t = np.arange(N) / fs
betas = [0.5, 1.0, 2.405, 5.0]

def fm_signal(beta):
    """s(t) = exp(j*beta*sin(2*pi*fm*t)): instantaneous frequency = beta*fm*cos(2*pi*fm*t)"""
    return np.exp(1j * beta * np.sin(2 * np.pi * fm * t))

def lines(s, nmax=12):
    """amplitude of the spectral line at n*fm, n = -nmax..nmax"""
    S = np.fft.fft(s) / N
    k = (np.arange(-nmax, nmax + 1) * int(fm)) % N
    return np.abs(S[k])

def bw98(a):
    """smallest +-n*fm band holding 98 % of the power"""
    p = a**2; c = len(a) // 2
    for n in range(c + 1):
        if p[c - n:c + n + 1].sum() >= 0.98 * p.sum():
            return 2 * n * fm
    return np.nan

def discriminator(s):
    """instantaneous frequency in Hz from the phase step between samples"""
    return np.angle(s[1:] * np.conj(s[:-1])) * fs / (2 * np.pi)

# 2. Measure every beta
print(f"{'beta':>6} {'df(kHz)':>8} {'J0 meas':>8} {'J0 th':>8} {'J1 meas':>8} {'J1 th':>8} "
      f"{'B98(kHz)':>9} {'Carson':>7} {'df meas':>8}")
rows = []
for b in betas:
    s = fm_signal(b); a = lines(s)
    c = len(a) // 2
    fi = discriminator(s)
    rows.append((b, a))
    print(f"{b:6.3f} {b*fm/1e3:8.3f} {a[c]:8.4f} {abs(jv(0, b)):8.4f} {a[c+1]:8.4f} {abs(jv(1, b)):8.4f} "
          f"{bw98(a)/1e3:9.1f} {2*(b*fm + fm)/1e3:7.1f} {fi.max()/1e3:8.3f}")

# 3. Carrier-level sweep: the carrier line follows |J0(beta)|
bsweep = np.arange(0, 6.01, 0.25)
j0m = np.array([lines(fm_signal(b), 0)[0] for b in bsweep])
err = np.max(np.abs(j0m - np.abs(jv(0, bsweep))))
print(f"Carrier sweep beta 0..6: max |J0 measured - J0 theory| = {err:.1e}")
print(f"First carrier null at beta = 2.405: carrier = {lines(fm_signal(2.405), 0)[0]:.4f}")

# 4. Figure
fig, ax = plt.subplots(2, 2, figsize=(11, 6.5))
for k, (b, col) in enumerate([(1.0, 'tab:blue'), (5.0, 'tab:purple')]):
    a = lines(fm_signal(b)); n = np.arange(-12, 13)
    ax[0, k].stem(n * fm / 1e3, a, linefmt=col, markerfmt='o', basefmt='k-')
    ax[0, k].axvspan(-(b * fm + fm) / 1e3, (b * fm + fm) / 1e3, color='orange', alpha=0.15, label="Carson's band")
    ax[0, k].set_title(f'Spectrum, beta = {b}'); ax[0, k].set_xlabel('Frequency offset (kHz)'); ax[0, k].set_ylabel('Line amplitude')
    ax[0, k].legend(loc='upper right'); ax[0, k].grid(alpha=0.3)
bb = np.linspace(0, 6, 400)
ax[1, 0].plot(bb, np.abs(jv(0, bb)), 'k-', label='|J0(beta)| theory')
ax[1, 0].plot(bsweep, j0m, 'o', color='tab:red', label='measured carrier')
ax[1, 0].set_xlabel('beta'); ax[1, 0].set_ylabel('Carrier amplitude'); ax[1, 0].set_title('Carrier line vs beta')
ax[1, 0].legend(); ax[1, 0].grid(alpha=0.3)
fi = discriminator(fm_signal(5.0))
tt = t[1:2001] * 1e3
ax[1, 1].plot(tt, fi[:2000] / 1e3, color='tab:green')
ax[1, 1].set_xlabel('Time (ms)'); ax[1, 1].set_ylabel('Instantaneous frequency (kHz)')
ax[1, 1].set_title('Discriminator output, beta = 5 (df = 5 kHz)'); ax[1, 1].grid(alpha=0.3)
fig.tight_layout(); fig.savefig('p3a_python.png', dpi=110)
# plt.show()
