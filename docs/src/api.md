# API Reference

## Mie series (Phase 0)

```@docs
mie_rcs_sphere
rcs_vs_size
optical_limit_rcs
wiscombe_nmax
```

## Model types

```@docs
AbstractRCSModel
MieSphere
rcs
rcs_vs_wavelength
```

## Bodies of revolution — closed-form RCS (Phase 1)

General PTD formula (paraboloid / spherical):

```@docs
scsn_general
```

Cone:

```@docs
cone_rcs
```

Canonical shape wrappers:

```@docs
paraboloid_rcs
spherical_segment_rcs
flat_disk_rcs
```
