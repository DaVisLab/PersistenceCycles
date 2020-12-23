# Persistence Cycles

<img src="./figures/pd_full.png" alt="drawing" height="120"/>
<img src="./figures/pairs_1perc.png" alt="drawing" height="120"/>
<img src="./figures/cycles.png" alt="drawing" height="120"/>
<img src="./figures/first_hurricane_cycles_zoom.png" alt="drawing" height="120"/>
<img src="./figures/full2-cycles-1perc.png" alt="drawing" height="120"/>




This repository contains a collection of modules for the Topology Toolkit (TTK) dedicated to computing persistent homology and representative cycles. A detailed description of the algorithms implemented in these modules can be found in the paper

.....




## Installation

You have two ways to install the modules provided in this repository. Both of them will require installing Paraview and TTK in your system. The [Topology ToolKit](https://topology-tool-kit.github.io) (TTK) is an open-source library and software collection for topological data analysis and visualization. 


### Installing TTK for the first time

If you have never used TTK before, it means you will have to compile Paraview and TTK from scratch following the instructions from the original [TTK website](https://topology-tool-kit.github.io/installation-0.9.7.html).

The only difference in the installation process is that you will have to use the TTK distribution provided in this repository (folder ttk-0.9.7) instead of downloading it from the TTK website.


### Integrating an existing TTK folder

If you already have TTK installed you can add the modules provided in this repository to your folder. You will have to create three modules and copy and past files contained in this repository.

Access yuor ttk folder and run the following commands

```
    ./scripts/createTTKmodule.sh BoundaryMatrix
    ./scripts/createTTKmodule.sh FormanGradient
    ./scripts/createTTKmodule.sh FG_PersistentHomology
```

The first line will create an empty module for implementing the [sequential algorithm](http://www.math.uchicago.edu/~shmuel/AAT-readings/Data%20Analysis%20/Edelsbrunner-Letscher-Zomordian.pdf) used for computing persistent homology. 
The second line will create an empty module for implementing the [Forman gradient](https://ieeexplore.ieee.org/document/5766002) at the base of the computation of persistent homology and the persistence cycles.
The third line will create an empty module for implementing the computation of the [persistence cycles]().

Once all modules have been created, you will have to remove files and folders that are not used in our plugin. You can do that by running the following command assuming you are in your TTK folder

```
    rm -R ./core/vtk/ttkBoundaryMatrix
    rm -R ./core/vtk/ttkFormanGradient
    rm -R ./paraview/BoundaryMatrix
    rm -R ./paraview/FormanGradient
    rm -R ./standalone/BoundaryMatrix
    rm -R ./standalone/FormanGradient
    rm -R ./standalone/FG_PersistentHomology
```

Now, you simply have to copy and past files provided in this repository to your local folder. Specifically, you have to substitute all files in the following folders

```
    ./core/base/boundaryMatrix
    ./core/base/formanGradient
    ./core/base/fG_PersistentHomology

    ./core/vtk/ttkFG_PersistentHomology

    ./paraview/FG_PersistentHomology
```

Once you are done, run the commands `cmake` and `make` again and the module `FG_PersistentHomology` should be part of your pool of modules in TTK.


## Visualizing persistence and the persistence pairs

Here we describe how to interact with the user interface in Paraview in order to visualize all the information computed by our module. Notice that you always have a chance to [run these modules automatically](https://topology-tool-kit.github.io/tutorials.html#python), for example in python, in order to automatize most of these steps.

### Output

![interface](./figures/outputs.png)

`FG_PersistentHomology` produces seven outputs. If you are not familiar with persistent homology, we suggest reading our [paper]() before continuing reading. The outputs produced by the plugin are:

- `Persistence Pairs` are the persistence pairs embedded in the domain of the input scalar field. Each pair is visualized with a line connecting a pair of points. The points correspond to the simplices creating and destroying a homology class. Each point is characterized by a `filtration value` and a `CelDimension` indicating the corresponding simplex dimension. Each line is characterized by a `filtration value` and a `type` indicating the homology class's dimension.

- `Persistence Diagram` the points populating the persistence diagram. Each point is characterized by a `filtration value` and a `Pair type` indicating the dimension of the homology class the point represents.

- `Homology` the points representing homology classes that never die. When working with images, this should contain a single point indicating the first vertex introduced in the filtration (i.e., global minimum).

- `1-cycles`, the persistence cycles computed for each persistence pair of type 1. Each persistence cycles is characterized by a unique identifier `CycleId` and a filtration value indicating the lifespan of such cycle

- `1-holes` persistence cycles are computed by visiting portions of the original dataset. In the case of persistence 1-cycles, these correspond to surfaces that are stored in `1-holes` for each persistence pair of type 1. Each hole is characterized by a unique identifier `CycleId` equal to the id of the corresponding cycle and a filtration value indicating such a cycle's lifespan.

- `2-cycles`, the persistence cycles computed for each persistence pair of type 2. Each persistence cycles is characterized by a unique identifier `CycleId` and a filtration value indicating the lifespan of such cycle

- `2-holes` persistence cycles are computed by visiting portions of the original dataset. In the case of persistence 2-cycles, these correspond to volumes that are stored in `2-holes` for each persistence pair of type 1. Each hole is characterized by a unique identifier `CycleId` equal to the id of the corresponding cycle and a filtration value indicating such a cycle's lifespan.



### Interface

![interface](./figures/interface.png)

Extraction of 1-cycles/2-cycles and 1-holes/2-holes is disabled by default and can be activated from the properties panel by clicking on the corresponding checkbox.

The user can also limit the extraction of persistence pairs, cycles, and holes only to homology classes with a certain lifespan. This is achieved by setting the interval of normalized persistence of interest (`Min Persistence`, `Max persistence`).