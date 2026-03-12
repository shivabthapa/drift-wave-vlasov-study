local Vlasov = G0.Vlasov
local Basis  = G0.Basis

pi = math.pi

-- -----------------------
-- Normalized constants
-- -----------------------
epsilon0 = 1.0
mu0      = 1.0

charge_elc = -1.0
charge_ion =  1.0
mass_elc   =  1.0
mass_ion   = 1.0   -- realistic-ish; you can reduce (e.g. 25, 100) for speed

-- Guide field (sets magnetization)
B0 = 1.0  -- along +z

-- -----------------------
-- Domain / resolution
-- -----------------------
Nx, Ny = 32, 32          -- start modest; increase after it runs
Nvx, Nvy, Nvz = 24, 24, 24

Lx, Ly = 20.0, 20.0      -- in your code units
vx_max, vy_max, vz_max = 6.0, 6.0, 6.0

poly_order = 2
basis_type = "serendipity"
time_stepper = "rk3"
cfl_frac = 0.9

t_end = 200.0
num_frames = 40

-- Seed mode (ky drives drift-wave-like structures)
ky_mode = 2.0*pi/Ly * 2.0
perturb_amp = 1.0e-6  -- small seed

-- -----------------------
-- Background equilibrium (periodic in x)
-- Use a double-tanh "bump" so n(x) matches at x=0 and x=Lx.
-- -----------------------
n0_bar = 1.0
delta_n = 0.3      -- strength of gradient/bump (keep <=0.3 to avoid violent transients)
w = 1.0            -- gradient width
x1 = 0.40*Lx
x2 = 0.60*Lx

-- Total pressure target (constant): p0 = n(x)*(Te(x)+Ti(x))
p0 = 1.0e-2        -- keep smallish initially to avoid huge thermal speeds
TeTi = 1.0         -- Te/Ti ratio (1 is fine to start)

local function sech2(z)
  local c = math.cosh(z)
  return 1.0/(c*c)
end

local function n0_profile(x)
  -- bump-shaped profile: edges match -> periodic OK
  local bump = 0.5*(math.tanh((x-x1)/w) - math.tanh((x-x2)/w))
  return n0_bar*(1.0 + delta_n*bump)
end

local function Te_profile(x)
  local n0 = n0_profile(x)
  -- pressure balance: Te+Ti = p0/n0
  local Tsum = p0/n0
  local Ti = Tsum/(1.0+TeTi)
  local Te = TeTi*Ti
  return Te
end

local function Ti_profile(x)
  local n0 = n0_profile(x)
  local Tsum = p0/n0
  local Ti = Tsum/(1.0+TeTi)
  return Ti
end

-- Maxwellian in 3V
local function maxwellian3v(n, m, T, vx, vy, vz, ux, uy, uz)
  local v2 = (vx-ux)*(vx-ux) + (vy-uy)*(vy-uy) + (vz-uz)*(vz-uz)
  local vt2 = T/m
  local norm = n / ((2.0*pi*vt2)^(1.5))
  return norm*math.exp(-v2/(2.0*vt2))
end

-- -----------------------
-- App
-- -----------------------
vlasovApp = Vlasov.App.new {
  tEnd = t_end,
  nFrame = num_frames,

  lower = {0.0, 0.0},
  upper = {Lx,  Ly},
  cells = {Nx,  Ny},

  cflFrac = cfl_frac,
  basis   = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- keep periodic as you requested
  periodicDirs = {1,2},
  decompCuts = {1},  -- fine for interactive tests

  -- -----------------------
  -- Electrons (now 3V)
  -- -----------------------
  elc = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_elc, mass = mass_elc,

    lower = {-vx_max, -vy_max, -vz_max},
    upper = { vx_max,  vy_max,  vz_max},
    cells = { Nvx,     Nvy,     Nvz},

    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.Func,
        init = function (t, xn)
          local x, y, vx, vy, vz = xn[1], xn[2], xn[3], xn[4], xn[5]
          local n0 = n0_profile(x)
          local Te = Te_profile(x)

          -- tiny seed perturbation (ky mode)
          local seed = 1.0 + perturb_amp*math.cos(ky_mode*y)

          -- Start with zero bulk flow (simplest). This is not a perfect equilibrium,
          -- but with small delta_n it should run and then you can refine.
          return seed * maxwellian3v(n0, mass_elc, Te, vx, vy, vz, 0.0, 0.0, 0.0)
        end
      }
    },

    evolve = true,
    diagnostics = {G0.Moment.M0, G0.Moment.M1}
  },

  -- -----------------------
  -- Ions (3V)
  -- -----------------------
  ion = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_ion, mass = mass_ion,

    lower = {-vx_max, -vy_max, -vz_max},
    upper = { vx_max,  vy_max,  vz_max},
    cells = { Nvx,     Nvy,     Nvz},

    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.Func,
        init = function (t, xn)
          local x, y, vx, vy, vz = xn[1], xn[2], xn[3], xn[4], xn[5]
          local n0 = n0_profile(x)
          local Ti = Ti_profile(x)

          -- same seed so quasi-neutral perturbation initially
          local seed = 1.0 + perturb_amp*math.cos(ky_mode*y)

          return seed * maxwellian3v(n0, mass_ion, Ti, vx, vy, vz, 0.0, 0.0, 0.0)
        end
      }
    },

    evolve = true,
    diagnostics = {G0.Moment.M0, G0.Moment.M1}
  },

  -- -----------------------
  -- EM Field (evolving)
  -- -----------------------
  field = Vlasov.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    init = function (t, xn)
      local x, y = xn[1], xn[2]

      -- Start with a uniform guide field Bz = B0 and no E.
      local Ex, Ey, Ez = 0.0, 0.0, 0.0
      local Bx, By, Bz = 0.0, 0.0, B0

      -- Last two are the error-cleaning fields (keep 0)
      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,

    evolve = true,
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

vlasovApp:run()
