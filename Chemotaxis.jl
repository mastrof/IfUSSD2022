using BacterialMotility
using LinearAlgebra
using Distributions
using Plots


#== system geometry ==#
L = 800 # box edge (μm)
d = 3 # dimensionality


#== convenience functions ==#
randposition() = sign.(rand(d) .- 1/2) .* (3L/8)

function randsphere(R)
    q = 2.0
    z = zeros(2)
    while q > 1.0
        x = rand(2)
        z .= 2 .* x .- 1
        q = z[1]*z[1] + z[2]*z[2]
    end # while
    s = sqrt(1.0 - q)
    [2z[1]*s, 2z[2]*s, 1.0-2.0*q] .* R
end # function


function randvelocity(d,U)
    x = rand(d) .- 1/2
    x ./ norm(x) .* U
end # function
randvelocity(U) = randvelocity(d,U)


#== custom motile patterns ==#
my_revflick!(b) = reverse_flick!(b, Normal(π,π/12), Uniform(3π/8, 5π/8))

#== population setup ==#
Δt = 0.1 # s
D_rot = 0.03 # rad²/s
s = propertiesBrumley("IntegrationTimestep" => Δt,
                      "RotationalDiffusivity" => D_rot,
                      "ReorientationRate" => 1.0,
                      "ChemotacticPrecision" => 0.5)

species = ["RRF"]
U = 40.0 # μm/s
N = 100
population = [BacteriumBrumley{d}(
    id = species[1], r = randsphere(0.9*L/2), v = randvelocity(U),
    turn! = my_revflick!, state = copy(s)) for _ in 1:N]

num_bacteria = length(population)


#== field setup ==#
R = 30.0 # μm
C = 50.0 # μM
Cbg = 0.01 # μM
f = Concentration_SteadyDiffusionSphericalSource_3D(R=R, C=C, Cbg=Cbg)

function save_traj!(traj, bs)
    nsteps, num_bacteria, d = size(traj)
    t = bs.clock[1]
    for i in 1:num_bacteria
        traj[t,i,:] .= bs.population[i].r
    end # for
end # function

nsteps = 700
trajectories = zeros(nsteps, num_bacteria, d)

function callback_outer(bs; kwargs...)
    save_traj!(trajectories, bs)
end # function


#== simulation ==#
bs = BacterialSystem(callback_outer = callback_outer,
                     population = population,
                     field = f)

integrate!(bs, nsteps)


#== plot setup ==#
plot_style(palette=:tab10) = (
    thickness_scaling = 1.5,
    guidefontsize = 12,
    tickfontsize = 12,
    legendfontsize = 8,
    grid = false,
    framestyle = :box,
    minorticks = true,
    tick_direction = :in,
    color_palette = palette,
    margin=3Plots.mm,
)

@userplot TailPlot
@recipe function f(p::TailPlot)
    traj, t1, t2, = p.args
    x = @view traj[t1:t2,:,1]
    y = @view traj[t1:t2,:,2]
    #z = @view traj[t1:t2,:,3]
    seriestype := :path
    linewidth --> 1
    label --> false
    x, y
end # recipe

@userplot HeadPlot
@recipe function f(p::HeadPlot)
    traj, t, = p.args
    x = @view traj[t:t,:,1]
    y = @view traj[t:t,:,2]
    #z = @view traj[t,:,3]
    seriestype := :scatter
    marker --> :circle
    markersize --> 5
    markerstrokewidth --> 0
    label --> false
    x, y
end # recipe

speciescolor = Dict(species .=> 1:length(species))
linecolor = [speciescolor[bact.id] for _ in 1:1, bact in population]

bgcolor = RGB(0.07, 0.07, 0.07)

function cfield(x,y)
    r = sqrt(x*x + y*y)
    r < R ? Cbg+C : Cbg+C*R/r
end # function
clims = (cfield(L/2,L/2), cfield(0,0))
xx = yy = -L/2:4:L/2
cc = cfield.(xx',yy)

#== plot final snapshot ==#
#=
colors = palette(:rainbow, num_bacteria)
q = plot(;plot_style()...)
α = range(0.25, 1.0; length=nsteps)
plot!(q, lims=(-L/2,L/2), axis=false, grid=false, size=(1200,1200))
heatmap!(q, xx, yy, cc, c=:bone, clims=clims, colorbar_ticks=false)
tailplot!(q, trajectories, 200, nsteps; lw=0.5, lc=(1:num_bacteria)')
headplot!(q, trajectories, nsteps; ms=8, mc=(1:num_bacteria)')
plot!(q, colorbar_title="Chemoattractant concentration",
      rightmargin=-14Plots.mm,
      colorbar_titlefontsize=20)
=#

#== plot animation ==#
ltail = 50
for t in 2:3:nsteps
    p = plot(;plot_style(:tab10)...)
    plot!(p, lims=(-L/2,L/2), aspect_ratio=1, axis=false, grid=false,
          bgcolor=bgcolor, size=(600,600))
    heatmap!(p, xx, yy, cc, c=:bone, clims=clims, colorbar=false)
    t0 = max(t-ltail, 1)
    tailplot!(p, trajectories, t0, t; lc=(1:num_bacteria)', lw=0.5)
    headplot!(p, trajectories, t; mc=(1:num_bacteria)', ms=4)
    #=
    for i in 1:length(species)
        plot!(p, [0.0], [0.0], lab=specieslabels[i],
              lc=speciescolor[species[i]],
              legend=:topright, legendfontsize=5,
              legendfontcolor=:white)
    end # for
    =#
#    δx = 100 # μm
#    xbar0 = -L/2 + L/30
#    xbar1 = xbar0 + δx
#    ybar = -L/2+L/40
#    plot!(p, [xbar0, xbar1], [ybar, ybar], lc=:white, lw=4, lab=false)
#    annotate!(p, xbar0+δx/2, ybar+22,
#              text("100 μm", :center, :white, 8))
#    xtime = L/2 - L/40
#    ytime = ybar
#    tnow = lpad(round(t*Δt, digits=1), 5, ' ')
#    annotate!(p, xtime, ytime,
#              text("t = $(tnow) s", :right, :white, 8))
    ndx = lpad(t-1, 4, '0')
    savefig(p, "frame$ndx.png")
end # for
