using Oceananigans
using Oceananigans.Units
using GLMakie

Lx = 15000kilometers
Ly = 9000kilometers

Nx = 75
Ny = 45

z_interfaces = [0.0,
                600.0,
                1300.0,
                2100.0,
                3000.0,
                4000.0,
                5100.0,
                6300.0,
                7600.0,
                8900.0,
                10300.0]

Nz = length(z_interfaces) - 1

grid = RectilinearGrid(size = (Nx, Ny, Nz),
                       x = (0, Lx),
                       y = (-Ly/2, Ly/2),
                       z = z_interfaces,
                       topology = (Periodic, Bounded, Bounded))

#model = NonhydrostaticModel(grid=grid)

coriolis = BetaPlane(f₀=1e-4, β=1.6e-11)
buoyancy = BuoyancyTracer()
advection = WENO()

#####
##### Build initial condition
#####

# Later, we'll want to introducing a forcing that relaxes
# The solution back to this initial condition, ie
#
# F = (bᵢ - b) / τ
#
# where τ is a parameter.

@inline ramp(y, Δy) = min(max(0, y/Δy + 1/2), 1)

N² = 4e-6 # [s⁻²] buoyancy frequency / stratification
M² = 8e-8 # [s⁻²] horizontal buoyancy gradient

Δy = Ly / 4       # width of the region of the front
Δb = Δy * M²      # buoyancy jump associated with the front
ϵb = 1e-2 * Δb    # noise amplitude

meridional_buoyancy(y) = Δb * ramp(y, Δy)
stratification(z) = N² * z
@inline bᵉ(x, y, z) = stratification(z) + meridional_buoyancy(y)
bᵢ(x, y, z) = bᵉ(x, y, z) + ϵb * randn()

model = NonhydrostaticModel(;
                            grid,
                            advection,
                            coriolis,
                            buoyancy,
                            tracers = :b,
                            timestepper = :RungeKutta3,
                           )

set!(model, b=bᵢ)

b = model.tracers.b

fig = Figure()
ax = Axis(fig[1, 1])
heatmap!(ax, interior(b, :, :, 4))

display(fig)

