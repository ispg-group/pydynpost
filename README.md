## What is pydynpost? 

pydynpost is a set of postprocessing scripts for both *ab initio* multiple spawning (AIMS) and trajectory surface hopping (TSH) methods.
The main motivation of this project is to gather all scripts necessary for data analysis of a previously performed nonadiabatic dynamics in one place. It features a high degree of modularity, allowing for inclusion of new functionality by simply adding a reference to it in the 
[aimsmodandmet.py](https://github.com/ispg-group/pydynpost/blob/master/aimsscripts/src/aimsmodandmet.py) / [tshmodandmet.py](https://github.com/ispg-group/pydynpost/blob/master/tshscripts/src/tshmodandmet.py) dictionary.

## What can it do?

The following table

| **Observables / metrics**            |      **AIMS**      |      **TSH**       |        **Task name** (*)          |
|-------------------------------------:|:------------------:|:------------------:|:---------------------------------:|
| State population                     | :heavy_check_mark: | :heavy_check_mark: |      `population`                 |
| Incoherent internals                 | :heavy_check_mark: | :heavy_check_mark: |      `internals`                  |
| Coherent internals                   | :heavy_check_mark: |     :question:     |      `internals`                  |
| Molecular population                 | :heavy_check_mark: | :heavy_check_mark: |      `molpop`                     |
| Comput. complexity & N<sub>TBF</sub> | :heavy_check_mark: |     :question:     |      `complexity`                 |
| Stochastic slections                 | :heavy_check_mark: |     :question:     |      `selection`                  |

:heavy_check_mark: = Accesible and implemented

:question: = Unaccesible or trivial

(*) This value is read into the `todo` variable. See section **Input parameters** for more information.

## Using pydynpost

pydynpost assumes that you are in a directory that contains the outputs of the nonadiabatic dynamics calculations you want to analyse.

   - For AIMS and TSH (without repetitions): `./${IC}m/`  (`${IC}` is arbitrary, e.g. `geom_`; `m` is an index for the initial condition)
   - For AIMSWISS / SSAIMS and TSH (with repetitons):  `./${RNG}n/${IC}m/`, (`${RNG}` is arbitrary, e.g. `rng`; `n` is an index for the run)
   
From within this working directory the postprocessing script can be launched by typing

```
  pydynpost
```

after which you enter interactive mode, in which the input parameters are parsed. This is handy for quick-and-dirty calculations, but inefficient when done more than once. For this reason pydynpost can also read the input parameters from the `dynpost.inp` file. This input file may look something like this:

```
dynMethod       = tsh
code            = ABIN
todo            = population
geomDir         = geom_
nrStates        = 6
nrRNGs          = 0
sampleSize      = 2
maxTime         = 188.19
step            = 0.1
```

## Input parameters

### General:

|  **Input variable name**    |       **Description**                |   **Possible options** | 
|:----------------------------|:-------------------------------------|:----------------------:|
| `dynMethod`                 |  Dynamics method                     |      `tsh` or `aims`   |
| `pckg` (only in aims)       |  Electronic structure interface      |      `tc` or `molpro`  |
| `code` (only in tsh)        |  Implementation of TSH               |        `ABIN`          |
| `todo`                      |  Observable / metric to be calcuated |      see (*)           |
| `sampleSize`                |  Number of intial conditions (ICs)   |       any `int`        |
| `nrRNGs`                    |  Number of runs per IC               |       any `int`        |


### Directory names: 

|  **Input variable name**           |       **Description**          | 
|:-----------------------------------|:-------------------------------|
| `outputDir` (only for aims+molpro) |   molpro subdirectory          |
| `geomDir`                          |   initial condition directory  |
| `RNGDir`                           |   rng seed directory           |

