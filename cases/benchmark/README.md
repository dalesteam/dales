DALES benchmark cases
=====================

RICO
----

Advection scheme 2 for mom, tke, thl, qt. 7 for sv
Bulk microphysics, no radiation.

### Modifications

* without custom surface scheme specified in RICO, instead using isurf=2
  (Forced surface temperature; fluxes are calculated)
* cloud field statistics off
* sampling statistics off
* cross section output on
* horizontal extent changed from 128^2 to 144^2, for flexible paralellizing options

### Variants

```
namoptions-144.001    300 seconds 1 core                  120.4 s on laptop
namoptions-fixed.001  100 seconds 1 core   fixed dt=2.5 s 133.4 s on laptop
namoptions-432.001    300
namoptions-864.001    300 
```

### Notes

No clouds or rain, time is too short

The fixed time step variant is useful for comparing older DALES versions,
the adaptive time stepping was made more strict in version 4.2.


RCEMIP
------

RRTMG radiation run every minute
simpleice(2) microphysics
advection scheme 5, except for qt: 77
504x504x146 cells, 200 m resolution.

This is the rce_small_les case with SST=300K from the RCE-MIP model intercomparison,
[Wing et al 2018](https://doi.org/10.5194/gmd-11-793-2018).

### Modifications

* fielddump off
* radstat off, because it assumes the radfull scheme and we are using RRTMG
* supposed to be run for a long time with constant solar angle, here started at noon: xtime=43200

### Variants

```
namoptions-96.300    96x96x146    300 s  1 core   92.7 s on laptop   1.3G RAM
namoptions-504.300   504x504x146  300 s                             ~35 G RAM
```

### Notes

RRTMG requires data files rrtmg_lw.nc, rrtmg_sw.nc. They are currently not distributed with DALES,
but can be obtained this way:

```
wget http://palm-model.org/trac/export/4419/palm/trunk/LIB/rrtmg/data/rrtmg_lw.nc
wget http://palm-model.org/trac/export/4419/palm/trunk/LIB/rrtmg/data/rrtmg_sw.nc

sha256sum rrtmg_*.nc
5ea5b5d4129d31952404985dc52059573f5d674cd8ee20bf5038f5e929f514a1  rrtmg_lw.nc
224aeba825ece5c53eebe97620a89fd37db6c85287aead2c624fccb80f36974e  rrtmg_sw.nc
```


Ruisdael
--------

RRTMG radiation, bulk (warm) microphysics, with large-scale forcings.
vertical grid stretched from 20 m to 150 m, up to 13.5 km


### Variants

```
namoptions-96.001 
96x96x160    15x15 km  156 m resolution  300 s   1 core  96.4 s on laptop     1.8 G RAM

namoptions-192.001
192x192x160   15x15 km  78 m resolution                                         5 G RAM

namoptions-1728.001,
1728x672x160 172x67 km 100 m resolution
```

1728 = 64*27
 672 = 32*3*7

### Notes

See above for obtaining RRTMG data.
Requires DALES 4.3 or above, due to a change in large scale forcing, ltimedepuv flag.


cao
---

A cold air outbreak case from Stephan de Roode et al.
Uses RRTMG radiation. Special for this case is that it
generates clouds and precipitation from the start, so
it can be used to measure how these impact the run time and
distribution over the subroutines.


GOAMAZON14
----------

A case with atmospheric chemistry and RRTMG radiation, from Jordi Vila et al.
Also included in cases/GOAMAZON14 with reference results.

### Variants

```
namoptions.047    36x36x200    300 s
```
