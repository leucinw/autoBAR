# HFE as a target in ForceBalance

Add target into forcebalance input file like this example:

```
$target
name HFE_CCl4
type Solvation_TINKER
weight 20.0
sfedata_txt hfe.txt
run_dynamic_every_iter 0
$end
```

Note:
`run_dynamic_every_iter`: 
- `0`: it will not do a full round FEP for iter_0001 and beyond to save time/
- `1`: it will do full round FEP at each iteration.

In `targets/` folder, add a folder like this example:

```
HFE_CCl4
├── gas.xyz
├── gas.key
├── hfe.txt
├── liquid.xyz
├── liquid.key
└── settings.yaml
```

In `hfe.txt`, simply put:
```
#HFE in kcal/mol
CCl4 0.08
```
