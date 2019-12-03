# cryoRayshader

The code enables users to visualize spatial data in the cryosphere in 3D animations, including a model of the glacier as well as graphical data visualization. The code is based on the original Rayshader (https://github.com/tylermorganwall/rayshader).
Necessary input data is a DEM of the area, spatial data of the glacier surface changing in time (e.g. data retrieved from satellites or models) and data to visualize on top including glacier velocities, surface features like lakes etc.

Note that distributed glacier thickness data exists for all glaciers globally (consult: https://www.glims.org/RGI/rgi60_dl.html). 
Velocity datasets also become more and more accesible (see for example https://nsidc.org/data/golive/)

Prerequisites are a working installation of R Studio (https://rstudio.com/) as well as working packages as detailed below.

We use data from our research studies to show examples of how the model can be used below.

Visualize modelled mass loss over a clean ice and a debris-covered glacier
------
Model outputs from [Wijngaard et al. (2019)](https://www.frontiersin.org/articles/10.3389/feart.2019.00143/full) (data available on request) allows us to visualize change in space. We do that below for a clean ice (Hintereis, European Alps) and debris-covered glacier (Langtang, Himalaya).

Simple mass loss data is visualized with orthophotos draped on top and the extent of the respective ice cover enhanced.

![](https://github.com/fidelsteiner/cryoRayshader/blob/master/exampleViz/double_example.gif)

Visualize modelled mass loss and associated data
------
In combination with the 3D visualization the data can be visualized alongside as a graph.

![](https://github.com/fidelsteiner/cryoRayshader/blob/master/exampleViz/massLossLangtangdata_example.gif)

Additionally spatial data like mass loss or velocities can be draped over the glacier surface.

![](https://github.com/fidelsteiner/cryoRayshader/blob/master/exampleViz/velocityLangtangdata_example.gif)

Visualize changing pond cover on a debris-covered glacier
------
We use pond outlines from [Steiner et al. (2019)](https://www.cambridge.org/core/journals/journal-of-glaciology/article/supraglacial-ice-cliffs-and-ponds-on-debriscovered-glaciers-spatiotemporal-distribution-and-characteristics/BEE84C3FF7F8BE25709171E8AE3BED5A) (data available https://doi.pangaea.de/10.1594/PANGAEA.899171) to show the change of water surfaces on a debris-covered glacier in the Himalaya.

Options include showing the numbers of features ...
![](https://github.com/fidelsteiner/cryoRayshader/blob/master/exampleViz/ponds_example.gif)

as well as the total and average area covered by the features ...

![](https://github.com/fidelsteiner/cryoRayshader/blob/master/exampleViz/ponds_example2.gif)
