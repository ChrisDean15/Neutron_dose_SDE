import numpy as np
import matplotlib.pyplot as plt
import itk

MeV_to_Gy = 1.602 * 10 ** (-13)


# Obtain SDE output from the txt file generated
def retrieve_sde_output(filename, airgap=0, bin_size=0.1):
    data = np.genfromtxt(filename)
    # Filter to region of interest, adjusted by airgap on x
    data = data[
        (data[:, 0] >= airgap)
        & (data[:, 0] < airgap + 18)
        & (data[:, 1] >= -5)
        & (data[:, 1] < 5)
        & (data[:, 2] >= -5)
        & (data[:, 2] < 5)
    ]
    x, y, z, vals = data.T
    x = x - airgap
    x_bins = int(18 / bin_size)
    y_bins = int(18 / bin_size)
    z_bins = int(18 / bin_size)
    volume = np.zeros((x_bins, y_bins, z_bins))

    for i in range(len(x)):
        xi = int(np.round(x[i] / bin_size))
        yi = int(np.round((y[i] + 5) / bin_size))  # shift y to [0,10]
        zi = int(np.round((z[i] + 5) / bin_size))  # shift z to [0,10]
        volume[xi, yi, zi] += vals[i]

    return volume


# Obtain geant4 output from the csv file generated
def retrieve_g4_output(filename):
    x, y, z, vals = np.genfromtxt(filename, unpack=True, delimiter=",", skip_header=1)
    x = x.astype(int)
    y = y.astype(int)
    z = z.astype(int)
    volume = np.zeros((100, 100, 100))
    for i in range(len(x)):
        volume[x[i], y[i], z[i]] = vals[i]

    return volume


# Compares central slice in a defined axis for both geant4 and SDE
def plot_slice(g4_mat, sde_mat, axis="x"):
    if axis == "x":
        g4_slice = g4_mat[49, :, :]
        sde_slice = sde_mat[49, :, :]
    elif axis == "y":
        g4_slice = g4_mat[:, 49, :]
        sde_slice = sde_mat[:, 49, :]
    else:
        g4_slice = g4_mat[:, :, 49]
        sde_slice = sde_mat[:, :, 49]
    # For a common colorbar ignoring infinite values when taking the log
    with np.errstate(divide="ignore"):
        g4_log = np.log10(g4_slice)
        sde_log = np.log10(sde_slice)

    g4_log[np.isneginf(g4_log)] = np.nan
    sde_log[np.isneginf(sde_log)] = np.nan

    vmin = np.nanmin([g4_log, sde_log])
    vmax = np.nanmax([g4_log, sde_log])

    f, ax = plt.subplots(1, 2, sharey=True, figsize=(11, 5))
    im1 = ax[0].imshow(
        g4_log.T, extent=[0, 10, -5, 5], origin="lower", vmin=vmin, vmax=vmax
    )
    im2 = ax[1].imshow(
        sde_log.T, extent=[0, 10, -5, 5], origin="lower", vmin=vmin, vmax=vmax
    )
    f.colorbar(
        im2,
        ax=ax.ravel().tolist(),
        label=r"Log(Dose[$\mu$Gy])",
        fraction=0.05,
        pad=0.04,
    )
    if axis == "x":
        ax[0].set_ylabel("y (cm)")
        ax[0].set_xlabel("z (cm)")
        ax[1].set_xlabel("z (cm)")
    elif axis == "y":
        ax[0].set_ylabel("z (cm)")
        ax[0].set_xlabel("x (cm)")
        ax[1].set_xlabel("x (cm)")
    else:
        ax[0].set_ylabel("y (cm)")
        ax[0].set_xlabel("x (cm)")
        ax[1].set_xlabel("x (cm)")
    ax[0].set_title("Geant4", fontweight="bold")
    ax[1].set_title("SDE", fontweight="bold")
    plt.suptitle(f"Central slice in {axis}", fontweight="bold")
    # plt.subplots_adjust(top = 1, bottom = 0.1, right = 0.9, left = 0.1, hspace = 0.3, wspace = 0.1)
    plt.show(block=True)


# Compares 2D projections
def plot_projection(g4_mat, sde_mat):
    g4_dose_proj = np.sum(g4_mat, axis=2)
    sde_dose_proj = np.sum(sde_mat, axis=2)

    # For a common colorbar ignoring infinite values when taking the log
    with np.errstate(divide="ignore"):
        g4_log = np.log10(g4_dose_proj)
        sde_log = np.log10(sde_dose_proj)

    g4_log[np.isneginf(g4_log)] = np.nan
    sde_log[np.isneginf(sde_log)] = np.nan

    vmin = np.nanmin([g4_log, sde_log])
    vmax = np.nanmax([g4_log, sde_log])

    f, ax = plt.subplots(1, 2, sharey=True, figsize=(11, 5))
    im1 = ax[0].imshow(
        g4_log.T, extent=[0, 10, -5, 5], origin="lower", vmin=vmin, vmax=vmax
    )
    im2 = ax[1].imshow(
        sde_log.T, extent=[0, 10, -5, 5], origin="lower", vmin=vmin, vmax=vmax
    )
    f.colorbar(
        im2,
        ax=ax.ravel().tolist(),
        label=r"Log(Dose[$\mu$Gy]), summed over z",
        fraction=0.05,
        pad=0.04,
    )
    ax[0].set_ylabel("y (cm)")
    ax[0].set_xlabel("x (cm)")
    ax[1].set_xlabel("x (cm)")
    ax[0].set_title("Geant4", fontweight="bold")
    ax[1].set_title("SDE", fontweight="bold")
    # plt.subplots_adjust(top = 1, bottom = 0.1, right = 0.9, left = 0.1, hspace = 0.3, wspace = 0.1)
    plt.show(block=True)


# Compares 1D dose either summed (projection) or from a slice
def plot_bragg_peak(
    g4_mat, sde_mat, how="proj"
):  # Define if projected Bragg peak or from a slice
    if how == "proj":
        g4_dose1d = np.sum(g4_mat, axis=(1, 2))
        sde_dose1d = np.sum(sde_mat, axis=(1, 2))
    else:
        sde_dose1d = sde_mat[:, 49, 49]
        g4_dose1d = g4_mat[:, 49, 49]
    print(sde_mat[101])
    x = np.linspace(0.1, 18, 180)
    print(len(x))
    print(len(sde_dose1d))
    plt.figure()
    plt.plot(x, sde_dose1d, label="SDE")
    plt.legend()
    plt.xlabel("Depth (cm)")
    plt.ylabel(r"Dose ($\mu$Gy)")
    plt.show(block=True)


# Shows all central slices for one 3D dose distribution only
def plot_all_slices(mat):
    slicex = mat[49, :, :]
    slicey = mat[:, 49, :]
    slicez = mat[:, :, 49]
    # For a common colorbar ignoring infinite values when taking the log
    with np.errstate(divide="ignore"):
        slicex_log = np.log10(slicex)
        slicey_log = np.log10(slicey)
        slicez_log = np.log10(slicez)

    slicex_log[np.isneginf(slicex_log)] = np.nan
    slicey_log[np.isneginf(slicey_log)] = np.nan
    slicez_log[np.isneginf(slicez_log)] = np.nan

    vmin = np.nanmin([slicex_log, slicey_log, slicez_log])
    vmax = np.nanmax([slicex_log, slicey_log, slicez_log])

    f, ax = plt.subplots(1, 3, sharey=True, figsize=(20, 7))
    im1 = ax[0].imshow(
        slicex_log.T, extent=[-5, 5, -5, 5], origin="lower", vmin=vmin, vmax=vmax
    )
    im2 = ax[1].imshow(
        slicey_log.T, extent=[0, 10, -5, 5], origin="lower", vmin=vmin, vmax=vmax
    )
    im3 = ax[2].imshow(
        slicez_log.T, extent=[0, 10, -5, 5], origin="lower", vmin=vmin, vmax=vmax
    )
    f.colorbar(
        im3,
        ax=ax.ravel().tolist(),
        label=r"Log(Dose[$\mu$Gy])",
        fraction=0.01,
        pad=0.04,
    )
    ax[0].set_title("X central slice", fontweight="bold")
    ax[1].set_title("Y central slice", fontweight="bold")
    ax[2].set_title("Z central slice", fontweight="bold")
    plt.subplots_adjust(
        top=1, bottom=0.1, right=0.88, left=0.05, hspace=0.05, wspace=0.1
    )
    plt.show(block=True)


# ---------------------------TODO: EXECUTION --------------------------------------
geant4 = retrieve_g4_output("Output/G4_TotalEneDepMap_1E6.csv")
g4_dose = geant4 * MeV_to_Gy *2.5* 10**5  # conversion to dose in microGy

sde = retrieve_sde_output("Output/test.txt")
sde_dose = sde # conversion to dose in microGy


plot_bragg_peak(g4_dose, sde_dose, how="nproj")
