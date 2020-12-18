# ...The Tea Is A Lie
### ...The Tea Is False
#### ...The Tea = 0
##### ...d Tea = 0
###### ...d t = 0
...steady state

This is a C++ software for calculating membrane concentration-profiles, fluxes, permeabilities, and so forth in steady-state based on Mesh Analysis. The software is built around a 2D rectangular membrane which is laterally periodic and therefore infinite, suitable for for example brick and mortar systems. The program is fundamentally constructed to model mass diffusion however it can also be used to model temperature flux and change diffusion. The major restriction of the scheme is that it assumes one constant concentration profile at the upper boundary of the system while utilizing a zink, i.e. zero concentration, at the lower boundary. For thermal/electric diffusion the temperature/electric potential is constant at those boundaries. Included is also the non-steady state software Advacuum which gives concentration-profiles and flows. That software is in many ways compatible with TTIAL, for example much of the input data is similar.

## Usage

### Requirements

- C++14 compiler
- The Eigen matrix library

### Building

The CMake build will automatically download Eigen.

~~~ bash
cmake .
make
~~~

### Examples

#### Steady State

A steady state example of a Brick and Mortar system can be run by executing the following commands

~~~ bash
cd examples/steadystate
../../ttial input.txt
~~~

Alternatively it is possible to use the Jupyter Notebook `brick_and_mortar.ipynb` in the same folder.

~~~ bash
jupyter-notebook brick_and_mortar.ipynb
~~~

Finally, in order to use IPython to run the software and produce figures simply convert the ipynb-file and run the produced py-file.

~~~ bash
jupyter nbconvert --to script brick_and_mortar.ipynb
ipython brick_and_mortar.py
~~~

#### Non Steady State

~~~ bash
cd examples/nonsteadystate
../../advacuum input.txt
~~~


### Input data

All input to the software is retrieved from the `input.txt` file which must be in the same folder as where the program is executed from. In the following we define the diffusion coefficient as D=U*R*T where U is the mobility, R the gas constant, and T the temperature. The solubility is defined by S = c0 / gamma / exp(mu0/RT) where c0 is the standard concentration, gamma the activity coefficient, and mu0 the standard chemical potential.
For charge diffusion `c_out` divited by `S_out` is the electric potential difference over the system and the resistances are the products between the solubilities `S` and diffusion coefficients `D`. 
For heat diffusion `c_out` divited by `S_out` is the temperature difference over the system and the thermal conductivities are the products between the solubilities `S` and diffusion coefficients `D`.

Important! Note that the input-file must end at the last input-row. If empty rows are present at the end of the file the program will crash.

Input parameter    |   Unit   |   Type   | Description
------------------ | -------- | -------- | -------------------
`model_nbr`	   | unitless | integer  | model number, 0 gives Brick and Mortar
`d`		   | m        | double   | width of bricks
`s`		   | m        | double   | horizontal spacing between bricks
`g`		   | m        | double   | vertical spacing between bricks
`t`		   | m        | double   | brick thickness
`N`		   | unitless | integer  | nbr of layers of bricks
`omega`		   | unitless | double   | offset ratio, negative gives random
`c_out`		   | kg/m^3   | double   | concentration outside system at upper boundary
`S_out`		   | kg/m^3   | double   | solubility outside system
`S_mv`		   | kg/m^3   | double   | mortar solubility (verticle)
`S_mh`		   | kg/m^3   | double   | mortar solubility (horizontal)
`S_bv`		   | kg/m^3   | double   | brick solubility (verticle)
`S_bh`		   | kg/m^3   | double   | brick solubility (horizontal)
`D_mv`		   | m^2/s    | double   | mortar diffusion coefficient (verticle)
`D_mh`		   | m^2/s    | double   | mortar diffusion coefficient (horizontal)
`D_bv`		   | m^2/s    | double   | brick diffusion coefficient (verticle)
`D_bh`		   | m^2/s    | double   | brick diffusion coefficient (horizontal)
`seed`		   | unitless | integer  | seed, negative gives random
`Nc`		   | unitless | integer  | number of columns of nodes
`Nr`		   | unitless | integer  | number of rows of nodes
`height`	   | m        | double   | height of system, only relevant when loading an external mesh
`width`		   | m        | double   | width of system, only relevant when loading an external mesh
`load_external`	   | N/A      | bool     | load external mesh, default: false
`input_folder`     | N/A      | string   | folder from which input is loaded, default: current directory, only relevant when loading an external mesh 
`output_folder`    | N/A      | string   | folder in which output is put, default:  current directory
`output_file`      | N/A      | string   | name of output-file, default: output.txt

#### Exclusive non-steady state input

Input parameter    |   Unit   |   Type   | Description
------------------ | -------- | -------- | -------------------
`time_steps`	   | unitless | integer  | number of iterations in time
`sample`	   | unitless | integer  | interval for output samples
`dt`		   | s        | double   | time-step
`time_periodic`	   | s        | double   | time between applying `c_out` at upper boundary
`evaporate`	   | N/A      | bool     | evaporate volume at upper boundary, default: false
`display`	   | N/A      | bool     | display progressbar or not, default: true

### Output data

Output files                |   Unit    | Description
--------------------------- | --------- | -------------
`output.txt`                |  -        | Miscellaneous data

Data in `output.txt`        |   Unit    | Description
--------------------------- | --------- | -------------
`model`                     |  unitless | Name of used model
`K_{M_ver/out}`             |  unitless | Partition coefficient between mortar (vertical) and outside
`K_{B_ver/out}`             |  unitless | Partition coefficient between bricks (vertical) and outside
`K_{M_ver/B_ver}`           |  unitless | Partition coefficient between mortar and bricks (both vertical)
`K_{M_hor/B_hor}`           |  unitless | Partition coefficient between mortar and bricks (both horizontal)
`height`                    |  m        | Height of system
`width`                     |  m        | Width of system
`Time`                      |   -       | Time to run the software [Days/Hours/Minutes/Seconds]

#### Exclusive steady state output

Output files                |   Unit    | Description
--------------------------- | --------- | -------------
`conc_hor_matrix.txt`       |  kg/m^3   | Concentration in the horizontal cells
`conc_ver_matrix.txt`       |  kg/m^3   | Concentration in the vertical cells
`sh_matrix.txt`             |  kg/m^3   | Solubility in the horizontal cells
`sv_matrix.txt`             |  kg/m^3   | Solubility in the vertical cells
`j_hor_matrix.txt`          |  kg/m^2 s | Flux in the horizontal cells
`j_ver_matrix.txt`          |  kg/m^2 s | Flux in the vertical cells
`j_hor_vector.txt`          |  kg/m^2 s | Mean of amplitude of flux in each row of the horizontal cells
`j_ver_vector.txt`          |  kg/m^2 s | Sum of flux in each row of the vertical cells
`Dh_matrix.txt`             |  m^2/s    | Diffusion coefficients in the horizontal cells
`Dv_matrix.txt`             |  m^2/s    | Diffusion coefficients in the vertical cells
`V_matrix.txt`              |  unitless | Potential in nodes in between cells

Data in `output.txt`        |   Unit    | Description
--------------------------- | --------- | -------------
`I_bot_sum`                 |  kg/m^2 s | Flux out of the system
`I_top_sum`                 |  kg/m^2 s | Flux into the system
`Reff`                      |  s m / kg | Effective resistance of the whole system
`j_ver`                     |  kg/m^2 s | Mean value of data in `jv_vector.txt`
`j_hor`                     |  kg/m^2 s | Mean value of data in `jh_vector.txt`

#### Exclusive non-steady state output

Output files                |   Unit    | Description
--------------------------- | --------- | -------------
`sn_matrix.txt`             |  kg/m^3   | Solubility in the nodes [FIX no such for steay-state]
`conc_X.txt`                |  kg/m^3   | Concentration profile, sample `X`
`volume_change.txt`         |  m^3      | volume of elements outside of upper boundary, vector in time
`mass_out.txt`              |  kg       | mass flowing out of the lower boundary, vector in time
