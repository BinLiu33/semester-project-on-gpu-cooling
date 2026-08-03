# Five-run comparison: tipsy mean profiles vs analytic (Sod).
# Example: file_keyword sod_10 -> sod_10_nocool, sod_10_kwh, sod_10_kwheat, sod_10_grackle, sod_10_grackleq
import os
import numpy as np
import pandas as pd
from sys import argv, exit
from matplotlib import pyplot as plt
import matplotlib.animation as animation
import sod

RUN_SUFFIXES = ("nocool", "kwh", "kwheat")


RUN_COLORS = {
    "nocool": "#1f77b4",    # blue
    "kwh": "#d62728",       # red
    "kwheat": "#e377c2",    # pink
    "grackleq": "#17becf",  # cyan
    "grackle": "#2ca02c",   # green
}

header_type = np.dtype(
    [('time', '=f8'), ('N', '=i4'), ('Dims', '=i4'), ('Ngas', '=i4'), ('Ndark', '=i4'), ('Nstar', '=i4'),
     ('pad', '=i4')])
gas_type = np.dtype(
    [('mass', '=f4'), ('x', '=f4'), ('y', '=f4'), ('z', '=f4'), ('vx', '=f4'), ('vy', '=f4'), ('vz', '=f4'),
     ('rho', '=f4'), ('temp', '=f4'), ('hsmooth', '=f4'), ('metals', '=f4'), ('phi', '=f4')])
dark_type = np.dtype(
    [('mass', '=f4'), ('x', '=f4'), ('y', '=f4'), ('z', '=f4'), ('vx', '=f4'), ('vy', '=f4'), ('vz', '=f4'),
     ('eps', '=f4'), ('phi', '=f4')])
star_type = np.dtype(
    [('mass', '=f4'), ('x', '=f4'), ('y', '=f4'), ('z', '=f4'), ('vx', '=f4'), ('vy', '=f4'), ('vz', '=f4'),
     ('metals', '=f4'), ('tform', '=f4'), ('eps', '=f4'), ('phi', '=f4')])

kB = 1.38064852e-23
u = 1.660538921e-27
M = 3.90e34
L = 6.17e17
T = 3e14
RHO = M / (L ** 3)
RHO *= 6.022e26 * 1e-6
L *= 3.240756e-17

if len(argv) < 6:
    print("Usage: python ./plot_compare5.py <pathName> <file_keyword> <step_interval> <nsteps> <temperature>")
    print("  file_keyword e.g. sod_10 -> loads sod_10_nocool, sod_10_kwh, ...")
    exit(1)

pathName = argv[1]
file_keyword = argv[2]
step_interval = int(argv[3])
nsteps = int(argv[4])
temp = int(argv[5])

if step_interval < 1:
    print("Error: step_interval must be >= 1")
    exit(1)

if temp == 10:
    temp = 9
    print("average =9, left = 10, right = 8")
if temp == 16000:
    temp = 14400
    print("average =14400, left = 16000, right = 12800")
if temp == 9:
    dt = 0.01
elif temp == 14400:
    dt = 0.00025
else:
    print("The temperature must be 9 or 14400 (or aliases 10 / 16000)")
    exit(1)

ach_names = [f"{file_keyword}_{suf}" for suf in RUN_SUFFIXES]

frame_steps = list(range(0, nsteps + 1, step_interval))
frame_cache = {}


def load_one_run(step_index, achOutName):
    key = (step_index, achOutName)
    if key in frame_cache:
        return frame_cache[key]

    file_step = step_index
    if file_keyword == "blast_part_2":
        file_step *= 10
    t = dt * file_step
    tipsy_path = os.path.join(pathName, achOutName, f"{achOutName}.{str(file_step).zfill(5)}")
    print("Loading %s" % tipsy_path)
    try:
        tipsy = open(tipsy_path, "rb")
    except OSError:
        print("Error: file not found %s" % tipsy_path)
        frame_cache[key] = None
        return None

    header = np.fromfile(tipsy, dtype=header_type, count=1)
    header = dict(zip(header_type.names, header[0]))
    gas = np.fromfile(tipsy, dtype=gas_type, count=header["Ngas"])
    _ = np.fromfile(tipsy, dtype=dark_type, count=header["Ndark"])
    _ = np.fromfile(tipsy, dtype=star_type, count=header["Nstar"])
    tipsy.close()

    x = gas["x"]
    rho = gas["rho"]
    temp_arr = gas["temp"]

    nbins = 1000
    bins = np.linspace(-0.5, 0.5, nbins + 1)
    count, _ = np.histogram(x, bins)
    rho_sum, _ = np.histogram(x, bins, weights=rho)
    temp_sum, _ = np.histogram(x, bins, weights=temp_arr)
    meanrho = np.divide(rho_sum, count, out=np.zeros_like(rho_sum), where=count > 0)
    meanT = np.divide(temp_sum, count, out=np.zeros_like(temp_sum), where=count > 0)
    meanPr = meanT * meanrho
    mean_x = (0.5 * (bins[:-1] + bins[1:])) * L

    time_years = t * T * 3.1689e-8
    out = {
        "file_step": file_step,
        "time_years": time_years,
        "mean_x": mean_x,
        "mean_density": meanrho * RHO,
        "mean_pressure": meanPr * kB / u * RHO,
        "mean_temperature": meanT,
    }
    frame_cache[key] = out
    return out


def analytic_at_step(step_index):
    file_step = step_index
    if file_keyword == "blast_part_2":
        file_step *= 10
    t = dt * file_step
    if "sod" not in file_keyword:
        return None
    (xgrid, PrE, _, rhoE, _, _, tempE) = sod.sod(t * T, temp)
    return {
        "x": xgrid[0, :],
        "density": rhoE[0, :],
        "pressure": PrE[0, :],
        "temperature": tempE[0, :],
    }


def preload_all():
    for frame_step in frame_steps:
        for ach in ach_names:
            load_one_run(frame_step, ach)


def render_mode(mode):
    fig = plt.figure()
    ax = plt.axes()
    colors = [RUN_COLORS[suf] for suf in RUN_SUFFIXES]

    def animate(step_index):
        ax.cla()
        frames = []
        for ach, c in zip(ach_names, colors):
            fr = load_one_run(step_index, ach)
            frames.append((ach, fr, c))
        if any(f[1] is None for f in frames):
            return

        first = next(f for f in frames if f[1] is not None)[1]
        for ach, fr, c in frames:
            if mode == "density":
                ax.plot(fr["mean_x"], fr["mean_density"], "-", linewidth=0.9, color=c, label=ach)
            elif mode == "pressure":
                ax.plot(fr["mean_x"], fr["mean_pressure"], "-", linewidth=0.9, color=c, label=ach)
            else:
                ax.plot(fr["mean_x"], fr["mean_temperature"], "-", linewidth=0.9, color=c, label=ach)

        an = analytic_at_step(step_index)
        if an is not None:
            if mode == "density":
                ax.plot(an["x"], an["density"], "k-", linewidth=0.75, label="analytic")
            elif mode == "pressure":
                ax.plot(an["x"], an["pressure"], "k-", linewidth=0.75, label="analytic")
            else:
                ax.plot(an["x"], an["temperature"], "k-", linewidth=0.75, label="analytic")

        if "sod" in file_keyword or "peq" in file_keyword:
            ax.set_xlim(-0.25 * L, 0.25 * L)
        else:
            ax.set_xlim(-0.5 * L, 0.5 * L)
        ax.set_xlabel(r"x [$pc$]")

        if mode == "density":
            ax.set_ylim(0.0, 2 * RHO)
            ax.set_ylabel(r"density [$AMU/cm^3$]")
        elif mode == "pressure":
            ax.set_ylim(0.0, 2 * temp * kB / u * RHO)
            ax.set_ylabel(r"pressure [$Pa$]")
        else:
            ax.set_ylim(0.0, 2 * temp)
            ax.set_ylabel(r"temperature [$K$]")

        ax.legend(loc="upper right", fontsize=7)
        ax.set_title(
            f"{file_keyword} (5 runs) STEP {first['file_step']:05d} TIME {first['time_years']:e} yr"
        )

    anim = animation.FuncAnimation(fig, func=animate, frames=frame_steps)
    gif_writer = animation.PillowWriter(fps=int(10 / step_interval))
    out_dir = pathName
    os.makedirs(out_dir, exist_ok=True)
    output_gif = os.path.join(out_dir, f"{file_keyword}_compare3_{mode}.gif")
    anim.save(output_gif, writer=gif_writer, dpi=150)
    plt.close(fig)


print("Preloading frames for all runs...")
preload_all()

for mode in ("density", "pressure", "temperature"):
    print(f"Rendering mode: {mode}")
    render_mode(mode)
