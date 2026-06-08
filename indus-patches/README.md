# INDUS patches

Patches to the INDUS PLUMED plugin used by this workflow.

## Applying

```bash
cd /path/to/indus-source
git am /path/to/ssi_workflow/indus-patches/*.patch
# rebuild PLUMED with the patched INDUS source
```

## 0001 — nearest union mode and configurable derivative sharing

**`union_mode = nearest`** for the union-of-spheres probe volume. The upstream
code computes the smooth indicator via the inclusion–exclusion product over all
overlapping spheres (`union_mode = product`). This patch adds a `nearest` mode
that computes the indicator from the distance to the nearest sphere center only.
For a hard-sphere union the membership criterion is identical; only the smooth
switching region (set by `sigma`, `alpha_c`) differs at boundaries. It avoids
the O(k) product over k overlapping spheres per water and is substantially
faster for large observation volumes.

```
ProbeVolume = {
  type        = union_spheres
  union_mode  = nearest      # or product (default)
  ...
}
```

**`NO_SHARE_ALL_DERIVATIVES`** flag for the order-parameters interface. When the
bias is applied through an external PLUMED `RESTRAINT` on `indus.ntilde`, INDUS
must broadcast all order-parameter derivatives to every MPI rank. That broadcast
is unnecessary when using the internal INDUS bias (`indus.ubias`). The flag lets
it be suppressed:

```
indus: INDUS INPUTFILE=... NO_SHARE_ALL_DERIVATIVES
```

Do **not** combine `NO_SHARE_ALL_DERIVATIVES` with a PLUMED `RESTRAINT` on
`indus.ntilde` — forces would be incorrect. This workflow uses the flag only in
the `internal-indus` bias mode.

---

**Note:** The included `.patch` is an abridged reference covering the key
interface change. Regenerate the authoritative patch from your INDUS working
tree with `git format-patch` before publishing, so line numbers and context
match the exact upstream commit you build against.
