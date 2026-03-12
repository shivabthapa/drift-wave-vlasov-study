-- Drift-wave instability test with the Hasegawa-Wakatani system.
-- Minimal modification of the standard turbulence example:
--   1) keep background density gradient drive
--   2) replace Gaussian blob with a small single-mode perturbation
--   3) shorten runtime to focus on linear / weakly nonlinear phase

local Vlasov = G0.Vlasov
local HasegawaWakatani = G0.Vlasov.Eq.HasegawaWakatani

alpha = 2.0   -- Adiabatic coupling constant
phi0  = 1.0e-4 -- Small seed amplitude

-- Single-mode perturbation
kx0 = 2.0*math.pi/40.0   -- one wavelength across x-box
ky0 = 2.0*math.pi/40.0   -- one wavelength across y-box

-- Simulation parameters.
Nx = 32
Ny = 32
Lx = 40.0
Ly = 40.0
poly_order = 2
basis_type = "serendipity"
time_stepper = "rk3"
cfl_frac = 1.0

-- Shorter run so we study instability before full turbulence dominates
t_end = 40.0
num_frames = 20

field_energy_calcs = GKYL_MAX_INT
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
    equation = HasegawaWakatani.new { alpha = alpha, is_modified = false },

    -- Background density gradient drive.
    -- Keeping the original linear profile is the least invasive choice.
    n0 = function(t, xn)
      local x, y = xn[1], xn[2]
      return x
    end,

    -- Small single-mode perturbation instead of Gaussian blob.
    -- This is better for tracking drift-wave growth and propagation.
    init = function(t, xn)
      local x, y = xn[1], xn[2]

      local phi  = phi0 * math.cos(kx0*x) * math.cos(ky0*y)
      local zeta = -(kx0*kx0 + ky0*ky0) * phi  -- zeta = ∇² phi

      return zeta, phi
    end,
  },

  skipField = true
}

vlasovApp:run()