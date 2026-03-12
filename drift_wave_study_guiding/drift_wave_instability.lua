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

-- Background magnetic field B0 z-hat.
B0 = 0.2

-- Background density parameters
n_bar    = 1.0
grad_amp = 0.1

-- Small electrostatic perturbation
phi0 = 5.0e-4
kx0  = 1.0
ky0  = 1.0

-- Domain
Lx = 2.0*pi
Ly = 2.0*pi

x0 = 0.5*Lx
kx_grad = 1.0

-- Grid
Nx, Ny   = 16, 16
Nvx, Nvy = 12, 12

vx_max = 3.0*vt_elc
vy_max = 3.0*vt_elc

poly_order   = 1
basis_type   = "serendipity"
time_stepper = "rk3"

cfl_frac = 0.4

t_end = 40.0
num_frames = 4

-- ------------------------------------------------------------
-- Density profile
-- ------------------------------------------------------------
local function n0_of_x(x)
  return n_bar * (1.0 + grad_amp * math.cos(kx_grad*(x - x0)))
end

local function dn0dx_of_x(x)
  return n_bar * (-grad_amp * kx_grad * math.sin(kx_grad*(x - x0)))
end

-- Electron pressure balance field
local function Ex0_of_x(x)
  local n0 = n0_of_x(x)
  return -(T_elc/n0) * dn0dx_of_x(x)
end

-- Perturbation potential
local function phi_of_xy(x,y)
  return phi0 * math.cos(kx0*x) * math.cos(ky0*y)
end

vlasovApp = Vlasov.App.new {

  tEnd = t_end,
  nFrame = num_frames,

  lower = {0.0,0.0},
  upper = {Lx,Ly},
  cells = {Nx,Ny},

  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  decompCuts = {1,1},
  periodicDirs = {1,2},

-- ============================================================
-- ELECTRONS
-- ============================================================

  elc = Vlasov.Species.new {

    modelID = G0.Model.Default,
    charge = charge_elc,
    mass   = mass_elc,

    lower = {-vx_max,-vy_max},
    upper = { vx_max, vy_max},
    cells = {Nvx,Nvy},

    numInit = 1,

    projections = {

      {
        projectionID = G0.Projection.Func,

        init = function(t,xn)

          local x,y,vx,vy = xn[1],xn[2],xn[3],xn[4]

          local Omega_e = charge_elc*B0/mass_elc

          -- guiding center
          local x_gc = x - vy/Omega_e

          -- periodic wrap
          x_gc = (x_gc % Lx)

          local n0 = n0_of_x(x_gc)

          local v2 = vx*vx + vy*vy

          local fM =
            (n0/(2*pi*vt_elc*vt_elc)) *
            math.exp(-v2/(2*vt_elc*vt_elc))

          return fM

        end
      }

    },

    evolve = true,

    diagnostics = {
      G0.Moment.M0,
      G0.Moment.M1
    }

  },

-- ============================================================
-- IONS
-- ============================================================

  ion = Vlasov.Species.new {

    modelID = G0.Model.Default,
    charge = charge_ion,
    mass   = mass_ion,

    lower = {-vx_max,-vy_max},
    upper = { vx_max, vy_max},
    cells = {Nvx,Nvy},

    numInit = 1,

    projections = {

      {
        projectionID = G0.Projection.Func,

        init = function(t,xn)

          local x,y,vx,vy = xn[1],xn[2],xn[3],xn[4]

          local Omega_i = charge_ion*B0/mass_ion

          local x_gc = x - vy/Omega_i

          x_gc = (x_gc % Lx)

          local n0 = n0_of_x(x_gc)

          local v2 = vx*vx + vy*vy

          local fM =
            (n0/(2*pi*vt_ion*vt_ion)) *
            math.exp(-v2/(2*vt_ion*vt_ion))

          return fM

        end
      }

    },

    evolve = true,

    diagnostics = {
      G0.Moment.M0,
      G0.Moment.M1
    }

  },

-- ============================================================
-- FIELD
-- ============================================================

  field = Vlasov.Field.new {

    epsilon0 = epsilon0,
    mu0      = mu0,

    init = function(t,xn)

      local x,y = xn[1],xn[2]

      local dphidx = -phi0*kx0*math.sin(kx0*x)*math.cos(ky0*y)
      local dphidy = -phi0*ky0*math.cos(kx0*x)*math.sin(ky0*y)

      local Ex_pert = -dphidx
      local Ey_pert = -dphidy

      local Ex = Ex0_of_x(x) + Ex_pert
      local Ey = Ey_pert
      local Ez = 0.0

      local Bx,By,Bz = 0.0,0.0,B0

      return Ex,Ey,Ez,Bx,By,Bz,0.0,0.0

    end,

    evolve = true,

    elcErrorSpeedFactor = 0.5,
    mgnErrorSpeedFactor = 0.5

  }

}

vlasovApp:run()