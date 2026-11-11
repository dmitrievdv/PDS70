using FITSIO
using GLMakie
using DataFrames
using Statistics
using Korg
using DataInterpolations

function rolling_sigma_clip(data, window, bound)
    clipped_data = deepcopy(data)
    median_data = deepcopy(data)
    # std_data = copy(data)
    n_data = length(data)
    for i = 1:n_data
        window_start = if i ≤ window
            1
        else
            i - window
        end

        window_end = if i ≥ n_data - window
            n_data
        else
            i + window
        end
        median_data[i] = median(data[window_start:window_end])
    end

    all_std =  std(data .- median_data, mean = 0.0)

    for i = 1:n_data
        # window_std = std(data[window_start:window_end] .- median_data[window_start:window_end], mean = 0.0)
        if abs(data[i] - median_data[i]) > bound*all_std
            println("Clipped! ", "$all_std $(median_data[i]) $(data[i])")
            clipped_data[i] = median_data[i]
        end 
    end
    return clipped_data
end

obs_dir = "archive_2" 

files = filter(name -> endswith(name, ".fits"), readdir("$obs_dir"))

fits_arr = FITS.("$obs_dir/" .* files)

v_rad = 6

dfs_Ha = map(fits_arr) do fits
    WAVE = vec(read(fits[2], "WAVE"))
    WAVE = WAVE / (1 + v_rad/3e5)
    FLUX = vec(read(fits[2], "FLUX"))
    DataFrame("WAVE" => WAVE[6538 .< WAVE .< 6588], "FLUX" => rolling_sigma_clip(FLUX[6538 .< WAVE .< 6588], 20, 3))
end

dfs_Hb = map(fits_arr) do fits
    WAVE = vec(read(fits[2], "WAVE"))
    WAVE .= WAVE / (1 + v_rad/3e5)
    FLUX = vec(read(fits[2], "FLUX"))
    DataFrame("WAVE" => WAVE[4836 .< WAVE .< 4886], "FLUX" => rolling_sigma_clip(FLUX[4836 .< WAVE .< 4886], 20, 3))
end



df_Ha_synth = begin
    fits=FITS("korg_ha.fits")

    WAVE = Korg.vacuum_to_air.(vec(read(fits[2], "WAVE")))
    
    FLUX = vec(read(fits[2], "FLUX"))
    CNTM = vec(read(fits[2], "CNTM"))

    # FLUX[6557 .< WAVE .< 6568] = CNTM[6557 .< WAVE .< 6568]

    interp_flux = AkimaInterpolation(FLUX, WAVE)
    interp_cntm = AkimaInterpolation(CNTM, WAVE)


    DataFrame("WAVE" => copy(dfs_Ha[1].WAVE), "FLUX" => interp_flux.(dfs_Ha[1].WAVE), "CNTM" => interp_cntm.(dfs_Ha[1].WAVE))
end

df_Hb_synth = begin
    fits=FITS("korg_hb.fits")
    WAVE = Korg.vacuum_to_air.(vec(read(fits[2], "WAVE")))
    FLUX = vec(read(fits[2], "FLUX"))
    CNTM = vec(read(fits[2], "CNTM"))

    interp_flux = AkimaInterpolation(FLUX, WAVE)
    interp_cntm = AkimaInterpolation(CNTM, WAVE)

    CNTM = vec(read(fits[2], "CNTM"))
    DataFrame("WAVE" => copy(dfs_Hb[1].WAVE), "FLUX" => interp_flux.(dfs_Hb[1].WAVE), "CNTM" => interp_cntm.(dfs_Hb[1].WAVE))
end

fig = Figure()
ax_Ha = Axis(fig[1,1], limits = ((6538,6588), nothing))

ax_Hb = Axis(fig[2,1], limits = ((4836,4886), nothing))

# deleteat!(dfs[5], 1)

for (i,df) in enumerate(dfs_Ha)
    lines!(ax_Ha, df.WAVE, df.FLUX ./ median(df.FLUX) / 1.01, label = files[i], alpha = 0.5)
end

for (i,df) in enumerate(dfs_Hb)
    lines!(ax_Hb, df.WAVE, df.FLUX ./ median(df.FLUX) / 1.2, label = files[i], alpha = 0.5)
end

Ha_phot_corr_prof = sum([df.FLUX / median(df.FLUX) / 1.01 for df in dfs_Ha])/length(files) - df_Ha_synth.FLUX ./ df_Ha_synth.CNTM .+ 1

lines!(ax_Ha, dfs_Ha[1].WAVE, Ha_phot_corr_prof, 
      label = "MEAN", color = :black, linewidth = 2)

renormed_synth_Ha = df_Ha_synth.FLUX ./ df_Ha_synth.CNTM
lines!(ax_Ha, df_Ha_synth.WAVE, renormed_synth_Ha)
vlines!(ax_Ha, [6562.81])

pruned_linelist_Ha = Korg.read_linelist("lines_Ha.h5")
pruned_wl_Ha = [Korg.vacuum_to_air(line.wl*1e8) for line in pruned_linelist_Ha]
pruned_sp_Ha = [line.species for line in pruned_linelist_Ha]
pruned_synth_ids_Ha = [findmin(x -> (x-wl)^2, df_Ha_synth.WAVE)[2] for wl in pruned_wl_Ha]
annotation!(ax_Ha, pruned_wl_Ha, fill(0.5, length(pruned_linelist_Ha)), 
                   pruned_wl_Ha, renormed_synth_Ha[pruned_synth_ids_Ha], 
                   text = string.(pruned_sp_Ha), labelspace = :data)


lines!(ax_Hb, dfs_Hb[1].WAVE, 
       sum([df.FLUX / median(df.FLUX) for df in dfs_Hb])/length(files) / 1.2 - df_Hb_synth.FLUX ./ df_Hb_synth.CNTM .+ 1,
       label = "MEAN", color = :black, linewidth = 2)
renormed_synth_Hb = df_Hb_synth.FLUX ./ df_Hb_synth.CNTM
lines!(ax_Hb, df_Hb_synth.WAVE, df_Hb_synth.FLUX ./ df_Hb_synth.CNTM )

pruned_linelist_Hb = Korg.read_linelist("lines_Hb.h5")
pruned_wl_Hb = [Korg.vacuum_to_air(line.wl*1e8) for line in pruned_linelist_Hb]
pruned_sp_Hb = [line.species for line in pruned_linelist_Hb]
pruned_synth_ids_Hb = [findmin(x -> (x-wl)^2, df_Hb_synth.WAVE)[2] for wl in pruned_wl_Hb]
annotation!(ax_Hb, pruned_wl_Hb[1:50], fill(0.3, length(pruned_linelist_Hb[1:50])), 
                   pruned_wl_Hb[1:50], renormed_synth_Hb[pruned_synth_ids_Hb[1:50]], 
                   text = string.(pruned_sp_Hb[1:50]), labelspace = :data)

axislegend(ax_Ha)

v_z_borders = [-6e7, 6e7]
λ_borders = 6562.81*(1 .+ v_z_borders/3e10)

Ha_prof_λs = dfs_Ha[1].WAVE[λ_borders[1] .< dfs_Ha[1].WAVE .< λ_borders[2]]
Ha_prof_vzs = 3e5*((Ha_prof_λs .- 6562.81)/6562.81)

Ha_prof_rs = Ha_phot_corr_prof[λ_borders[1] .< dfs_Ha[1].WAVE .< λ_borders[2]]

# lines(Ha_prof_vzs, Ha_prof_rs)

open("pds70obs-$obs_dir-nolines.dat", "w") do io
    for (v_z, r) in zip(Ha_prof_vzs, Ha_prof_rs)
        println(io, "$v_z $r")
    end
end

fig
