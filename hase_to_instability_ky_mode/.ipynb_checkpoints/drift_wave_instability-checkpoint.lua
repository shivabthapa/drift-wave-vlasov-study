-- Drift-wave instability test with the Hasegawa-Wakatani system.
-- Minimal changes from the standard example:
--   1) keep background density-gradient drive
--   2) use a small single-mode ky perturbation
--   3) initialize density-like perturbation out of phase with phi
--   4) run long enough to detect linear growth

local Vlasov = G0.Vlasov
local HasegawaWakatani = G0.Vlasov.Eq.HasegawaWakatani

alpha = 10.0      -- adiabatic coupling constant
phi0  = 1.0e-4   -- small seed amplitude

-- Use a pure ky mode (better aligned with drift-wave drive)
kx0 = 0.0
ky0 = 2.0*math.pi/10.0   -- one wavelength across y-box

-- Simulation parameters
Nx = 64
Ny = 64
Lx = 40.0
Ly = 40.0

poly_order   = 2
basis_type   = "serendipity"
time_stepper = "rk3"
cfl_frac     = 1.0

-- Run longer to capture instability growth
t_end = 120.0
num_frames = 20

field_energy_calcs   = GKYL_MAX_INT
integrated_mom_calcs = GKYL_MAX_INT
integrated_L2_f_calcs = GKYL_MAX_INT
dt_failure_tol = 1.0e-4
num_failures_max = 20

vlasovApp = Vlasov.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedL2fCalcs = integrated_L2_f_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,

  lower = { -Lx/2.0, -Ly/2.0 },
  upper = {  Lx/2.0,  Ly/2.0 },
  cells = { Nx, Ny },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  decompCuts = { 1, 1 },
  periodicDirs = { 1, 2 },

  fluid = Vlasov.FluidSpecies.new {
    equation = HasegawaWakatani.new {
      alpha = alpha,
      is_modified = false
    },

    -- Background density-gradient drive
    n0 = function(t, xn)
      local x = xn[1]
      return x
    end,

    -- Small drift-wave-like perturbation
    init = function(t, xn)
      local x, y = xn[1], xn[2]

      -- electrostatic potential perturbation
      local phi = phi0 * math.cos(ky0*y)

      -- vorticity zeta = ∇^2 phi
      local zeta = -(ky0*ky0) * phi

      -- second evolved field initialized out of phase
      -- this helps avoid starting in a nearly trivial aligned state
      local npert = phi0 * math.sin(ky0*y)

      return zeta, npert
    end,
  },

  skipField = true
}

vlasovApp:run()