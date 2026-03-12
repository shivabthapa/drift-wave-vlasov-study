-- Drift-wave turbulence test with the Hasegawa–Wakatani system

local Vlasov = G0.Vlasov
local HasegawaWakatani = G0.Vlasov.Eq.HasegawaWakatani

alpha = 10.0

-- larger amplitude to reach nonlinear regime
phi0  = 1.0e-2

-- grid resolution (slightly higher for turbulence)
Nx = 64
Ny = 64

Lx = 40.0
Ly = 40.0

poly_order   = 2
basis_type   = "serendipity"
time_stepper = "rk3"
cfl_frac     = 1.0

-- run longer to allow nonlinear saturation
t_end = 400.0
num_frames = 40

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

    -- background density gradient
    n0 = function(t, xn)
      local x = xn[1]
      return x
    end,

    -- broadband perturbation seed
    init = function(t, xn)

      local x, y = xn[1], xn[2]

      -- random potential perturbation
      local phi = phi0 * (math.random() - 0.5)

      -- start vorticity near zero
      local zeta = 0.0

      -- random density perturbation
      local npert = phi0 * (math.random() - 0.5)

      return zeta, npert
    end,
  },

  skipField = true
}

vlasovApp:run()