## What is pydynpost? 

pydynpost is a set of postproccessing scripts for both *ab initio* multiple spawning (AIMS) and trajectory surface hopping (TSH) methods.
The main motivation of this project is to gather all scripts necessary for data analysis of a previously performed nonadiabatic dynamics in one place. It features a high degree of modularity, allowing for inclusion of new functionality by simply adding a reference to it in the 
[aimsmodandmet.py](https://github.com/ispg-group/pydynpost/blob/master/aimsscripts/src/aimsmodandmet.py) / [tshmodandmet.py](https://github.com/ispg-group/pydynpost/blob/master/tshscripts/src/tshmodandmet.py) dictionary.

## What can it do?


| **Observables / metrics**            |      **AIMS**      |      **TSH**       |
|-------------------------------------:|:------------------:|:------------------:|
| State population                     | :heavy_check_mark: | :heavy_check_mark: |
| Incoherent internals                 | :heavy_check_mark: | :heavy_check_mark: |
| Coherent internals                   | :heavy_check_mark: |     :question:     |
| Molecular population                 | :heavy_check_mark: | :heavy_check_mark: |
| Comput. complexity & N<sub>TBF</sub> | :heavy_check_mark: |     :question:     |

:heavy_check_mark: = Accesible and implemented

:question: = Unaccesible or trivial

## Using pydynpost

pydynpost assumes that you are in a directory that contains the outputs of the nonadiabatic dynamics calculations you want to analyse.

   - For AIMS and TSH (without repetitions): `./${IC}m/`  (`${IC}` is arbitrary, e.g. `geom_`; `m` is an index for the initial condition)
   - For AIMSWISS / SSAIMS and TSH (with repetitons):  `./${RNG}n/${IC}m/`, (`${RNG}` is arbitrary, e.g. `rng`; `n` is an index for the run)
   
From within this working directory the postproccessing script can be launched by typing

```
  pydynpost
```

after which you enter interactive mode, in which the input parameters are parsed. This is handy for quick-and-dirty calculations, but inefficient when done more than once. For this reason pydynpost can also read the input parameters from the `dynpost.inp` file. This input file may look something like this:

```
dynMethod       = tsh
code            = ABIN
todo            = population
internalName    = Cr1-C2
dissPartners    = Cr-C
thresh          = 0.10335104202
internalType    = bl
geomDir         = geom_
nrStates        = 6
nrRNGs          = 0
sampleSize      = 2
maxTime         = 188.19
step            = 0.1
duplicatesExist = n
dupList         =
```

## Input parameters

