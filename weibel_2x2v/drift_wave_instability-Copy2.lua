local Vlasov = G0.Vlasov
pi = math.pi

-- Normalized code units.
epsilon0 = 1.0
mu0      = 1.0

-- Species parameters.
charge_elc = -1.0
mass_elc   = 1.0

charge_ion =  1.0
mass_ion   = 25.0

-- Temperatures (thermal speeds).
T_elc = 0.02
T_ion = 0.02

vt_elc = math.sqrt(T_elc/mass_elc)
vt_ion = math.sqrt(T_ion/mass_ion)

-- Background magnetic field B0 z-hat. (reduced to ease stiffness)
B0 = 0.2

-- Density gradient profile n0(x) = n_bar * [1 + a * tanh((x-x0)/L_n)].
n_bar = 1.0
grad_amp = 0.1

-- Small electrostatic potential perturbation:
phi0 = 5.0e-4
kx0  = 1.0
ky0  = 1.0

-- Domain
Lx = 2.0*pi
Ly = 2.0*pi

-- Gradient geometry
x0 = 0.5*Lx
Ln = 0.15*Lx

-- Very cheap grid
Nx, Ny   = 8, 8
Nvx, Nvy = 8, 8

-- Tighten velocity bounds (cheaper & less stiff)
vx_max = 3.0*vt_elc
vy_max = 3.0*vt_elc

poly_order   = 1
basis_type   = "serendipity"
time_stepper = "rk3"

-- Relax CFL to avoid dt-collapse cascades
cfl_frac = 0.5

t_end = 20.0
num_frames = 1

-- Background density
local function n0_of_x(x)
  return n_bar * (1.0 + grad_amp * math.tanh((x - x0)/Ln))
end

-- Seed potential
local function phi_of_xy(x,y)
  return phi0 * math.cos(kx0*x) * math.cos(ky0*y)
end

vlasovApp = Vlasov.App.new {
  tEnd = t_end,
  nFrame = num_frames,

  lower = { 0.0, 0.0 },
  upper = { Lx, Ly },
  cells = { Nx, Ny },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  decompCuts = { 1 },
  periodicDirs = { 1, 2 },

  -- -----------------
  -- Electrons
  -- -----------------
  elc = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_elc, mass = mass_elc,

    lower = { -vx_max, -vy_max },
    upper = {  vx_max,  vy_max },
    cells = { Nvx, Nvy },

    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.Func,
        init = function(t, xn)
          local x, y, vx, vy = xn[1], xn[2], xn[3], xn[4]
          local n0 = n0_of_x(x)
          local v2 = vx*vx + vy*vy
          local fM = (n0/(2*pi*vt_elc*vt_elc)) * math.exp(-v2/(2*vt_elc*vt_elc))
          return fM
        end
      }
    },

    evolve = true,
    diagnostics = { G0.Moment.M0, G0.Moment.M1 }
  },

  -- -----------------
  -- Ions (frozen for fast sanity run)
  -- -----------------
  ion = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_ion, mass = mass_ion,

    lower = { -vx_max, -vy_max },
    upper = {  vx_max,  vy_max },
    cells = { Nvx, Nvy },

    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.Func,
        init = function(t, xn)
          local x, y, vx, vy = xn[1], xn[2], xn[3], xn[4]
          local n0 = n0_of_x(x)
          local v2 = vx*vx + vy*vy
          local fM = (n0/(2*pi*vt_ion*vt_ion)) * math.exp(-v2/(2*vt_ion*vt_ion))
          return fM
        end
      }
    },

    evolve = false
  },

  -- -----------------
  -- Maxwell field
  -- -----------------
  field = Vlasov.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    init = function(t, xn)
      local x, y = xn[1], xn[2]

      -- Seed electrostatic E = -grad(phi).
      local dphidx = -phi0*kx0*math.sin(kx0*x)*math.cos(ky0*y)
      local dphidy = -phi0*ky0*math.cos(kx0*x)*math.sin(ky0*y)

      local Ex = -dphidx
      local Ey = -dphidy
      local Ez = 0.0

      local Bx, By, Bz = 0.0, 0.0, B0

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,

    evolve = true,

    -- Key change: enable error-speed factors to reduce stiffness / dt-collapse.
    elcErrorSpeedFactor = 1.0,
    mgnErrorSpeedFactor = 1.0
  }
}

vlasovApp:run()
