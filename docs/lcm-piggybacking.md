# Warm-cloud LCM piggybacking

This opt-in mode keeps warm bulk microphysics (`imicro=2`) as the driver of
winds, pressure, buoyancy and radiation. LCM follows the same flow while
maintaining separate prognostic vapor and potential-temperature tracers.
LCM phase changes heat/dry only that passenger state, never the bulk state.

```fortran
&nammicrophysics
 imicro = 2
 l_rain = .false.
 l_sedc = .false.
/
&namlcm
 lcm_piggyback = .true.
 switch_koehler_micro = .true.
 switch_analytical_micro = .false.
 switch_lem_micro = .false.
 switch_ice = .false.
 switch_collection = .false.
 switch_sedimentation = .false.
 switch_particle_source = .false.
 switch_sgs_velocities = .false.
/
```

Use the DALES executable with a double-precision CPU build, LCM enabled,
periodic cold starts and `ladaptive=.false.`. Open boundaries, restart initialization, independent
passenger radiative/surface schemes and nudging are not supported by this
first mode. The `daleslib` time loop is not wired for this mode. The existing
LCM-only selection `imicro=12` is retained.

The user-facing switch is `lcm_piggyback` in DALES's `NAMLCM`, defaulting to
false. It selects how the host couples LCM, so it is not duplicated in
`lcm_config_t`. Standalone LCM provides the same generic readouts and optional
refresh callback to any host; it does not know DALES's bulk scheme or own the
host's atmospheric fields.

The `external/lcm` submodule pins the matching feedback-interface commit.
Alternatively, build that LCM revision separately and configure DALES with
`-DENABLE_LCM=ON -DLCM_ROOT=/path/to/lcm/build`.

## Passenger fields

- `lcm_theta`: dry potential temperature, K. DALES transports this without
  applying its saturation adjustment. The adapter derives absolute
  temperature with the current Exner profile and passes dry static energy
  `T + g*z/cp` in LCM's temperature-unit convention.
- `lcm_qv`: vapor mass per unit DALES reference-air mass. Reference density
  `rhobf` matches the conservative scalar transport weights. The adapter
  converts this to vapor per dry-air mass for LCM using the passenger's own
  vapor plus represented liquid, not the bulk driver's diagnosed cloud water.
- `lcm_t_prev`, `lcm_p_prev`, `lcm_qv_prev`: transported previous passenger
  state used by the normal LCM non-LEM condensation calculation. These are
  refreshed after passenger phase-change feedback, including RK levels and
  halos.

The passenger receives the driver's surface heat/vapor fluxes and radiative
heating as shared external forcing. Prescribed horizontal/temporal sources
are shared. Scalar subsidence uses each passenger tracer's own gradient;
particles also receive the prescribed vertical subsidence velocity. Thermal
sponge damping uses each passenger field's own horizontal mean, with the
driver's height-dependent rate. It acts on `theta` and vapor, whereas the
driver damps `theta_l` and total water: these are equivalent variables only
in cloud-free air. This is a passenger sponge rule, not an exact copy of the
bulk operator in cloudy sponge cells. The three history tracers receive no
extra sponge tendency; their means are refreshed for DALES's top boundary.
Subgrid scalar diffusion is DALES's existing tracer diffusion. Particle SGS random
motion is disabled in this first configuration.

The two DALES RK arrays can have independently randomized cold-start values.
Both passenger levels are initialized from their corresponding DALES levels.
Initial particle liquid is included when matching total water and the
DALES constant-cp/constant-L moist-enthalpy invariant. Passenger fields are
not overwritten by bulk temperature or vapor on subsequent steps.

## Phase-change exchange

LCM refreshes attached passenger thermodynamics after particle migration and
sources, just before condensation. Its public mass readouts return kg/cell,
with particle multiplicity already included. The net condensation increment
excludes transported liquid; evaporation has a negative sign.

For cell reference-air mass `M = rhobf*dx*dy*dz` and net particle-water change
`dm`, the adapter applies only to passenger tracers:

```
qv_passenger -= dm/M
theta_passenger += L_DALES/(cp_DALES*Exner) * dm/M
```

This conserves `M*qv + particle_liquid` and the phase-change invariant
`cp*Exner*theta + L*qv`. It is not a demonstration that transport and all
external forcings conserve global moist energy by themselves.

Host phase feedback uses DALES's `L=2.53e6 J/kg`, `cp=1004 J/(kg K)` so the
bulk/passenger comparison has the same energy convention. LCM's unchanged
growth coefficients internally use `L=2.50e6 J/kg`; LCM retains its own
Magnus-type saturation law, the non-LEM supersaturation/relaxation sequence,
Koehler effects and original numerical Euler solver. No direct-S override or
safeguarded solver is introduced. These remaining scheme differences can
produce legitimate cloud differences. Equality of cloud liquid is not a
pass criterion.

The necessary `l_rain` setting is broadcast to every MPI rank. Constant
cold-start runs now honor `ladaptive=false` from their first full timestep,
including bulk-only controls; the adaptive startup ramp remains available
when adaptive stepping is enabled.

The Pa/hPa correction in LCM's growth coefficients is an independent physics
prerequisite, not a piggybacking algorithm change. It should be reviewed
separately from the host feedback interface. The rain broadcast and fixed-step
startup correction are DALES prerequisites shared with bulk-only controls.

Experimental diagnostics and run scripts are kept outside the repository.
