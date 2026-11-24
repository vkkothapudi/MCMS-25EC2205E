"""P3B: Phase modulation - PM vs FM, the narrowband phase modulator and PM demodulation (MCMS 25EC2205E) - Python version.
Compares the frequency deviation of PM and FM as the tone frequency changes, shows what a
triangle message does to each, measures the error of the narrowband (Armstrong) phase modulator
and recovers the message from the phase.  Needs: numpy, matplotlib."""
import numpy as np
import matplotlib.pyplot as plt

# 1. Parameters
fs = 200e3                       # sample rate (S/s)
N = int(fs)                      # 1 s of samples
t = np.arange(N) / fs
kp = 2.0                         # PM phase sensitivity (rad per volt)
kf = 5e3                         # FM frequency sensitivity (Hz per volt)
Am = 1.0                         # message amplitude (V)

def tone(f): return Am * np.cos(2 * np.pi * f * t)
def tri(f):  return Am * (2 * np.abs(2 * ((f * t + 0.25) % 1) - 1) - 1)        # triangle, +-Am
def pm(m):   return np.exp(1j * kp * m)                                         # phase follows m(t)
def fm(m):   return np.exp(1j * 2 * np.pi * kf * np.cumsum(m) / fs)             # frequency follows m(t)
def inst_freq(s): return np.angle(s[1:] * np.conj(s[:-1])) * fs / (2 * np.pi)  # Hz

# 2. Peak deviation vs tone frequency: FM fixed, PM grows with fm
print('Part 1: tone message, Am = 1 V, kp = 2 rad/V, kf = 5 kHz/V')
print(f"{'fm(Hz)':>7} {'dfFM meas':>10} {'dfFM th':>8} {'dfPM meas':>10} {'dfPM th':>8} {'betaFM':>7} {'betaPM':>7}")
fms = [250, 500, 1000, 2000, 4000]
dpm, dfm = [], []
for f in fms:
    m = tone(f)
    a, b = inst_freq(fm(m)).max(), inst_freq(pm(m)).max()
    dfm.append(a); dpm.append(b)
    print(f"{f:7d} {a:10.1f} {kf*Am:8.1f} {b:10.1f} {kp*Am*f:8.1f} {kf*Am/f:7.2f} {kp*Am:7.2f}")

# 3. Triangle message: PM frequency jumps between two levels, FM sweeps smoothly
ftri = 500
fp = inst_freq(pm(tri(ftri))); ff = inst_freq(fm(tri(ftri)))
th = kp * 4 * Am * ftri / (2 * np.pi)                # |dm/dt| = 4*Am*f for a triangle
mid = np.abs(fp) > 0.5 * th                          # samples away from the jumps
print(f"\nPart 2: triangle message, {ftri} Hz")
print(f"PM: frequency jumps between +-{np.median(np.abs(fp[mid])):.1f} Hz (theory +-{th:.1f} Hz), "
      f"{np.mean(np.abs(np.abs(fp) - th) < 0.01 * th)*100:.1f} % of samples within 1 % of a level")
print(f"FM: frequency sweeps smoothly between {ff.min():.0f} and {ff.max():.0f} Hz (theory +-{kf*Am:.0f} Hz)")

# 4. Narrowband phase modulator: cos(wc t) - theta*sin(wc t) at baseband = 1 + j*theta
print('\nPart 3: narrowband phase modulator, s = 1 + j*theta (theta = k*m, tone 1 kHz)')
print(f"{'theta pk':>8} {'phase pk':>9} {'error %':>8} {'th err %':>9} {'env max':>8} {'th env':>7}")
nb_k = [0.1, 0.2, 0.5, 1.0]
nb_err = []
for k in nb_k:
    theta = k * tone(1e3)
    s = 1 + 1j * theta
    ph = np.angle(s).max()
    e = (k - ph) / k * 100; nb_err.append(e)
    print(f"{k:8.2f} {ph:9.4f} {e:8.2f} {(k-np.arctan(k))/k*100:9.2f} {np.abs(s).max():8.4f} {np.sqrt(1+k*k):7.4f}")

# 5. PM demodulator: the phase of s(t) divided by kp gives m(t) back
sine = Am * np.sin(2 * np.pi * 1e3 * t)               # a sine starts at 0, so the phase starts unwrapped
m = sine * 1.4                                        # peak phase 2.8 rad, still inside +-pi
mh = np.angle(pm(m)) / kp
print(f"\nPart 4: PM demodulation, peak phase {kp*1.4:.1f} rad: rms error = {np.sqrt(np.mean((mh-m)**2)):.1e} V")
m2 = sine * 2.0                                       # peak phase 4 rad: needs unwrapping
mh2 = np.unwrap(np.angle(pm(m2))) / kp
print(f"Peak phase {kp*2.0:.1f} rad (beyond pi): wrapped error = {np.sqrt(np.mean((np.angle(pm(m2))/kp-m2)**2)):.2f} V, "
      f"after unwrap = {np.sqrt(np.mean((mh2-m2)**2)):.1e} V")

# 6. Figure
fig, ax = plt.subplots(2, 2, figsize=(11, 6.5))
n = slice(0, 800); tm = t[1:801] * 1e3
ax[0, 0].plot(tm, ff[n] / 1e3, color='tab:red', label='FM of triangle')
ax[0, 0].plot(tm, fp[n] / 1e3, color='tab:purple', label='PM of triangle')
ax[0, 0].set_xlabel('Time (ms)'); ax[0, 0].set_ylabel('Instantaneous frequency (kHz)')
ax[0, 0].set_title('Triangle message, 500 Hz'); ax[0, 0].legend(loc='upper right'); ax[0, 0].grid(alpha=0.3)
ax[0, 1].plot(fms, np.array(dfm) / 1e3, 'o-', color='tab:red', label='FM (measured)')
ax[0, 1].plot(fms, np.array(dpm) / 1e3, 's-', color='tab:purple', label='PM (measured)')
ax[0, 1].plot(fms, kp * Am * np.array(fms) / 1e3, 'k--', lw=1, label='PM theory kp*Am*fm')
ax[0, 1].set_xlabel('Tone frequency fm (Hz)'); ax[0, 1].set_ylabel('Peak deviation (kHz)')
ax[0, 1].set_title('Deviation vs tone frequency'); ax[0, 1].legend(); ax[0, 1].grid(alpha=0.3)
kk = np.linspace(0.02, 1.2, 100)
ax[1, 0].plot(kk, (kk - np.arctan(kk)) / kk * 100, 'k-', label='theory (theta - atan theta)/theta')
ax[1, 0].plot(nb_k, nb_err, 'o', color='tab:green', label='measured')
ax[1, 0].set_xlabel('Peak phase theta (rad)'); ax[1, 0].set_ylabel('Phase error (%)')
ax[1, 0].set_title('Narrowband phase modulator error'); ax[1, 0].legend(); ax[1, 0].grid(alpha=0.3)
tt = t[:600] * 1e3
ax[1, 1].plot(tt, m2[:600], 'k--', lw=1.2, label='message sent')
ax[1, 1].plot(tt, np.angle(pm(m2))[:600] / kp, color='tab:orange', lw=1, label='phase / kp (wrapped)')
ax[1, 1].plot(tt, mh2[:600], color='tab:blue', lw=1.6, alpha=0.8, label='unwrapped phase / kp')
ax[1, 1].set_xlabel('Time (ms)'); ax[1, 1].set_ylabel('Message (V)')
ax[1, 1].set_title('PM demodulation, peak phase 4 rad'); ax[1, 1].legend(loc='lower right', fontsize=8); ax[1, 1].grid(alpha=0.3)
fig.tight_layout(); fig.savefig('p3b_python.png', dpi=110)
# plt.show()
