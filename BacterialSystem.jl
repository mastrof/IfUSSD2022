using BacterialMotility
using LinearAlgebra
using Distributions
using Plots


#== system geometry ==#
L = 1.5e3 # box edge (μm)
d = 2 # dimensionality


#== convenience functions ==#
randposition() = (rand(d) .- 1/2) .* (L/2)
function randvelocity(d,U)
    x = rand(d) .- 1/2
    x ./ norm(x) .* U
end # function
randvelocity(U) = randvelocity(d,U)


#== custom motile patterns ==#
my_helical!(b) = tumble!(b, Normal(π/4, π/32))
my_tumble!(b) = tumble!(b, Uniform(-π, π))
my_revflick!(b) = reverse_flick!(b, Normal(π,π/12), Uniform(3π/8, 5π/8))
my_reverse!(b) = tumble!(b, Normal(π,π/16))

#== population setup ==#
Δt = 0.01 # s
D_rot = 0.03 # rad²/s
s1 = properties("IntegrationTimestep" => Δt,
                "RotationalDiffusivity" => D_rot,
                "ReorientationRate" => 1.3,)
s2 = properties("IntegrationTimestep" => Δt,
                "RotationalDiffusivity" => D_rot)
s3 = properties("IntegrationTimestep" => Δt,
                "RotationalDiffusivity" => D_rot,
                "ReorientationRate" => 0.75)
s4 = properties("IntegrationTimestep" => Δt,
                "RotationalDiffusivity" => D_rot)

species = ["RT_1", "RT_2", "RRF_1", "RR"]
U = 40.0 # μm/s
N = 1
r01 = [-1.0, 1.0] .* L/6
r02 = [1.0, 1.0] .* L/6
r03 = [-1.0, -1.0] .* L/6
r04 = [1.0, -1.0] .* L/6
pop1 = [Bacterium{d}(
    id = species[1], r = r01, v = randvelocity(U),
    run! = run!, turn! = my_helical!, state = copy(s1)) for _ in 1:N]
pop2 = [Bacterium{d}(
    id = species[2], r = r02, v = randvelocity(U),
    run! = run!, turn! = my_tumble!, state = copy(s2)) for _ in 1:N]
pop3 = [Bacterium{d}(
    id = species[3], r = r03, v = randvelocity(1.1*U),
    run! = run!, turn! = my_revflick!, state = copy(s3)) for _ in 1:N]
pop4 = [Bacterium{d}(
    id = species[4], r = r04, v = randvelocity(1.3*U),
    run! = run!, turn! = my_reverse!, state = copy(s4)) for _ in 1:N]

population = vcat(pop2)#, pop2, pop3, pop4)
num_bacteria = length(population)


function callback_inner(b,f; kwargs...)
    rotational_diffusion!(b)
    #pbc!(b)
end # function

function save_traj!(traj, bs)
    nsteps, num_bacteria, d = size(traj)
    t = bs.clock[1]
    for i in 1:num_bacteria
        traj[t,i,:] .= bs.population[i].r
    end # for
end # function

nsteps = 800
trajectories = zeros(nsteps, num_bacteria, d)

function callback_outer(bs; kwargs...)
    save_traj!(trajectories, bs)
end # function


#== simulation ==#
bs = BacterialSystem(callback_inner = callback_inner,
                     callback_outer = callback_outer,
                     population = population)

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
    x = @view traj[t,:,1]
    y = @view traj[t,:,2]
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

#== plot final snapshot ==#
q = plot(;plot_style(:Dark2)...)
plot!(q, axis=false, grid=false, size=(600,600))
tailplot!(q, trajectories, 1, nsteps; lc=1, lw=3)
headplot!(q, trajectories, nsteps; mc=1, ms=15)
display(q)

#== plot animation ==#
ltail = 100
for t in 2:4:0nsteps
    p = plot(;plot_style(:Dark2)...)
    plot!(p, lims=(-L/2,L/2), aspect_ratio=1, axis=false, grid=false,
          bgcolor=bgcolor, size=(600,600))
    t0 = max(t-ltail, 1)
    tailplot!(p, trajectories, t0, t; lc=linecolor)
    headplot!(p, trajectories, t; mc=permutedims(linecolor))
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
