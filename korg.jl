using Korg
using GLMakie
using DataFrames
using FITSIO

linelist = Korg.read_linelist("gfallvac08oct17.dat", format = "kurucz_vac")

linelist[1051558] = Korg.Line(6.571028e-5, -0.45, Korg.Species("Fe I"), 4.733140788520537)

A_X = format_A_X(0.0, Dict("Fe" => 0.1, "Ti" => -0.2, "Ca" => -0.1, "Al" => -2, "Co" => -2, "Cr" => -2))
atm = Korg.interpolate_marcs(4040, 4.44, A_X)
result = synthesize(atm, linelist, A_X, [(6500, 6600, 0.01)], use_MHD_for_hydrogen_lines=false)

λs = result.wavelengths
rot_flux = Korg.apply_rotation(result.flux, λs, 17)
continuum = result.cntm

data = Dict("WAVE" => λs,
        "FLUX" => rot_flux,
        "CNTM" => continuum)

FITS("korg_ha.fits", "w") do f
    write(f, data; hdutype = TableHDU)
end

pruned_linelist = Korg.prune_linelist(atm, linelist, A_X, (6535,6590,0.01), threshold = 0.1, sort_by_EW = true)
Korg.save_linelist("lines_Ha.h5", pruned_linelist)

# result = synthesize(atm, linelist, A_X, [(4800, 4900, 0.01)], use_MHD_for_hydrogen_lines=false)

# λs = result.wavelengths
# rot_flux = Korg.apply_rotation(result.flux, λs, 17)
# continuum = result.cntm

# data = Dict("WAVE" => λs,
#         "FLUX" => rot_flux,
#         "CNTM" => continuum)

# FITS("korg_hb.fits", "w") do f
#     write(f, data; hdutype = TableHDU)
# end

# pruned_linelist = Korg.prune_linelist(atm, linelist, A_X, (4840,4890,0.01), threshold = 0.1, sort_by_EW = true)
# Korg.save_linelist("lines_Hb.h5", pruned_linelist)
