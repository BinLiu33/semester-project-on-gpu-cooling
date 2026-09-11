# temperature of peq is 14400
# python plot_compare.py ./results sod_16000_nocool 1 100 14400
# plot pressure, density, and temperature gif of given file
import numpy as np
import pandas as pd
from sys import argv, exit
from matplotlib import pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import matplotlib.animation as animation
import sod

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

skip = 1
nsteps = int(argv[4])

kB = 1.38064852e-23  # Boltzmann constant in SI units
u = 1.660538921e-27  # atomic mass unit in kg
M = 3.90e34  # code unit mass in kg
L = 6.17e17  # code unit length in m
T = 3e14  # code unit time in s (T = L^3/(G*M))
RHO = M / (L ** 3)  # code unit density in kg/m^3
# convert from kg/m^3 to amu/m^3
RHO *= 6.022e26 * 1e-6
L *= 3.240756e-17  # convert from m to pc
temp = int(argv[5])
if temp == 10:
    temp = 9
    print("average =9, left = 10, right = 8")
if temp == 16000:
    temp = 14400
    print("average =14400, left = 16000, right = 12800")
if (temp == 9) :
    dt = 0.01
elif temp == 14400:
    dt = 0.00025 
else:
    print("The temperature must be 9 or 14400")
    exit(1)   

pathName = argv[1]

achOutName = argv[2]
step_interval = int(argv[3])

if len(argv) < 6:
    print("Usage: python ./plot_compare.py <pathName> <file_keyword> <step_interval> <nsteps> <temperature>")
    exit(1)
if step_interval < 1:
    print("Error: step_interval must be >= 1")
    exit(1)

mean = True
SCATTER_STRIDE = 10
GIF_DPI = 110
SCATTER_CMAP = "viridis"  
frame_steps = list(range(0, nsteps + 1, step_interval))
frame_cache = {}

def load_frame_data(step_index):
    if step_index in frame_cache:
        return frame_cache[step_index]

    file_step = step_index
    if achOutName == "blast_part_2":
        file_step *= 10

    t = dt * file_step

    filename = "%s/%s/%s.%s" % (
        pathName,
        achOutName,
        achOutName,
        str(file_step).zfill(5)
    )

    print("Loading", filename)

    try:
        tipsy = open(filename, "rb")
    except:
        print("Error: file not found")
        frame_cache[step_index] = None
        return None

    header = np.fromfile(
        tipsy,
        dtype=header_type,
        count=1
    )

    header = dict(zip(header_type.names, header[0]))

    gas = np.fromfile(
        tipsy,
        dtype=gas_type,
        count=header['Ngas']
    )

    _ = np.fromfile(
        tipsy,
        dtype=dark_type,
        count=header['Ndark']
    )

    _ = np.fromfile(
        tipsy,
        dtype=star_type,
        count=header['Nstar']
    )

    tipsy.close()

    # --------------------------
    # numpy arrays
    # --------------------------

    x = gas['x']
    rho = gas['rho']
    temp_arr = gas['temp']

    # --------------------------
    # sorting only for plotting
    # --------------------------

    order = np.argsort(rho)

    x_sorted = x[order]
    rho_sorted = rho[order]
    temp_sorted = temp_arr[order]

    
    x_plot = x_sorted[::SCATTER_STRIDE]
    rho_plot = rho_sorted[::SCATTER_STRIDE]
    temp_plot = temp_sorted[::SCATTER_STRIDE]
    color_plot = np.log10(rho_plot)  

    # --------------------------
    # fast binning
    # --------------------------

    nbins = 1000

    bins = np.linspace(
        -0.5,
        0.5,
        nbins + 1
    )

    count, _ = np.histogram(
        x,
        bins
    )

    rho_sum, _ = np.histogram(
        x,
        bins,
        weights=rho
    )

    temp_sum, _ = np.histogram(
        x,
        bins,
        weights=temp_arr
    )

    meanrho = np.divide(
        rho_sum,
        count,
        out=np.zeros_like(rho_sum),
        where=count > 0
    )

    meanT = np.divide(
        temp_sum,
        count,
        out=np.zeros_like(temp_sum),
        where=count > 0
    )

    meanPr = meanT * meanrho

    mean_x = (
        0.5 * (bins[:-1] + bins[1:])
    ) * L

    # --------------------------
    # analytic solution
    # --------------------------

    analytic = None

    if "sod" in achOutName:
        (
            xgrid,
            PrE,
            _,
            rhoE,
            _,
            _,
            tempE
        ) = sod.sod(
            t * T,
            temp
        )

        analytic = {
            "x": xgrid[0, :],
            "density": rhoE[0, :],
            "pressure": PrE[0, :],
            "temperature": tempE[0, :]
        }

    time_years = (
        t * T * 3.1689e-8
    )

    frame_cache[step_index] = {
        "file_step": file_step,
        "time_years": time_years,
        "x": x_plot * L,
        "rho": rho_plot,
        "temp": temp_plot,
        "color": color_plot,
        "mean_x": mean_x,
        "mean_density": meanrho * RHO,
        "mean_pressure": meanPr * kB / u * RHO,
        "mean_temperature": meanT,
        "analytic": analytic
    }

    return frame_cache[step_index]

def render_mode(mode):
    fig = plt.figure()
    ax = plt.axes()

    
    sm = ScalarMappable(norm=Normalize(vmin=COLOR_VMIN, vmax=COLOR_VMAX), cmap=SCATTER_CMAP)
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label(r'particle $\log_{10}\rho$ (code units)')

    def animate(step_index):
        frame = load_frame_data(step_index)
        if frame is None:
            return

        ax.cla()
        if mode == "density":
            ax.scatter(frame["x"], frame["rho"] * RHO, s=0.5, edgecolors="none", c=frame["color"],
                       cmap=SCATTER_CMAP, vmin=COLOR_VMIN, vmax=COLOR_VMAX)
        elif mode == "pressure":
            ax.scatter(frame["x"], frame["temp"] * frame["rho"] * kB / u * RHO, s=0.5, edgecolors="none",
                       c=frame["color"], cmap=SCATTER_CMAP, vmin=COLOR_VMIN, vmax=COLOR_VMAX)
        else:
            ax.scatter(frame["x"], frame["rho"] * RHO, s=0.5, edgecolors="none", c=frame["color"],
                       cmap=SCATTER_CMAP, vmin=COLOR_VMIN, vmax=COLOR_VMAX)

        if mean:
            if mode == "density":
                ax.plot(frame["mean_x"], frame["mean_density"], 'r-',
                        linewidth=0.75, label=r'mean density',antialiased=False)
            elif mode == "pressure":
                ax.plot(frame["mean_x"], frame["mean_pressure"], 'r-',
                        linewidth=0.75, label=r'mean pressure',antialiased=False)
            elif mode == "temperature":
                ax.plot(frame["mean_x"], frame["mean_temperature"], 'r-', linewidth=0.75,
                        label=r'mean temperature',antialiased=False)
            else:
                ax.plot(frame["mean_x"], frame["mean_density"], 'r-',
                        linewidth=0.75, label=r'mean density')


        if frame["analytic"] is not None:
            if mode == "density":
                ax.plot(frame["analytic"]["x"], frame["analytic"]["density"], 'k-', linewidth=0.75,
                        label=r'analytic density')
            elif mode == "pressure":
                ax.plot(frame["analytic"]["x"], frame["analytic"]["pressure"], 'k-', linewidth=0.75,
                        label=r'analytic pressure')
            elif mode == "temperature":
                ax.plot(frame["analytic"]["x"], frame["analytic"]["temperature"], 'k-', linewidth=0.75,
                        label=r'analytic temperature')


        if "sod" in achOutName:
            ax.set_xlim(-0.25 * L, 0.25 * L)
        elif "peq" in achOutName:
            ax.set_xlim(-0.25 * L, 0.25 * L)
        else:
            ax.set_xlim(-0.5 * L, 0.5 * L)
        ax.set_xlabel(r'x [$pc$]')

        if mode == "density":
            ax.set_ylim(0.0, 2 * RHO)
            ax.set_ylabel(r'density [$AMU/cm^3$]')
        elif mode == "pressure":
            ax.set_ylim(0.0, 2.5 * temp * kB / u * RHO)
            ax.set_ylabel(r'pressure [$Pa$]')
        elif mode == "temperature":
            ax.set_ylim(0.0, 2 * temp)
            ax.set_ylabel(r'temperature [$K$]')
        else:
            ax.set_ylim(0.0, 2 * RHO)
            ax.set_ylabel(r'density [$kg/m^3$]')

        ax.legend(loc='upper right')
        ax.set_title(f"{achOutName} STEP {frame['file_step']:05d} TIME {frame['time_years']:e}yr")

    anim = animation.FuncAnimation(fig, func=animate, frames=frame_steps)
    gif_writer = animation.PillowWriter(fps=int(10/step_interval))
    output_gif = f"{pathName}/{achOutName}/{achOutName}_compare_{mode}.gif"
    anim.save(output_gif, writer=gif_writer, dpi=GIF_DPI)
    plt.close(fig)

print("Preloading frames once for all modes...")
for frame_step in frame_steps:
    load_frame_data(frame_step)


_all_colors = np.concatenate(
    [f["color"] for f in frame_cache.values() if f is not None and len(f["color"]) > 0]
)
COLOR_VMIN = float(np.nanmin(_all_colors))
COLOR_VMAX = float(np.nanmax(_all_colors))

for mode in ("density", "pressure", "temperature"):
    print(f"Rendering mode: {mode}")
    render_mode(mode)

