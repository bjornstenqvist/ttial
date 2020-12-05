# ...The Tea Is A Lie
### ...The Tea Is False
#### ...The Tea = 0
##### ...d Tea = 0
###### ...d t = 0
...steady state

This is a C++ software for calculating membrane concentration-profiles, fluxes, permeabilities, and so forth in steady-state based on Mesh Analysis. The software is built around a 2D rectangular membrane which is laterally periodic and therefore infinite, suitable for for example brick and mortar systems. The program is fundamentally constructed to model mass diffusion however it can also be used to model temperature flux and change diffusion. The major restriction of the scheme is that it assumes one constant concentration profile at the upper boundary of the system while utilizing a zink, i.e. zero concentration, at the lower boundary. For thermal/electric diffusion the temperature/electric potential is constant at those boundaries.

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

### Example

~~~ bash
cd example
../ttial
~~~

### Input data

In the following we define the diffusion coefficient as D=U*R*T where U is the mobility, R the gas constant, and T the temperature. The solubility is defined by S = c0 / gamma / exp(mu0/RT) where c0 is the standard concentration, gamma the activity coefficient, and mu0 the standard chemical potential. 
For charge diffusion `c_out` divited by `S_out` is the electric potential difference over the system and the resistances are the products between the solubilities `S` and diffusion coefficients `D`. 
For heat diffusion `c_out` divited by `S_out` is the temperature difference over the system and the thermal conductivities are the products between the solubilities `S` and diffusion coefficients `D`.

Input parameter    |   Unit   |   Type   | Description
------------------ | -------- | -------- | -------------------
`model_nbr`	   | unitless | integer  | model number, 0 gives Brick and Mortar
`d`		   | m        | double   | width of bricks
`s`		   | m        | double   | horizontal spacing between bricks
`g`		   | m        | double   | vertical spacing between bricks
`t`		   | m        | double   | brick thickness
`N`		   | unitless | integer  | nbr of layers of bricks
`omega`		   | unitless | double   | offset ratio, negative gives random
`c_out`		   | kg/m^3   | double   | concentration outside system
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
`height`	   | m        | double   | height of system, only relevant when loading an external mesh
`width`		   | m        | double   | width of system, only relevant when loading an external mesh
`load_external`	   | N/A      | bool     | load external mesh, true/false
`Nc`		   | unitless | integer  | number of columns of nodes
`Nr`		   | unitless | integer  | number of rows of nodes

### Output data

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
`output.txt`                |  N/A      | Miscellaneous data

Data in `output.txt`        |   Unit    | Description
--------------------------- | --------- | -------------
`model`                     |  unitless | Name of used model
`K_{M_ver/out}`             |  unitless | Partition coefficient between mortar (vertical) and outside
`K_{B_ver/out}`             |  unitless | Partition coefficient between bricks (vertical) and outside
`K_{M_ver/B_ver}`           |  unitless | Partition coefficient between mortar and bricks (both vertical) 
`K_{M_hor/B_hor}`           |  unitless | Partition coefficient between mortar and bricks (both horizontal) 
`I_bot_sum`                 |  kg/m^2 s | Flux out of the system
`I_top_sum`                 |  kg/m^2 s | Flux into the system
`Reff`                      |  s m / kg | Effective resistance of the whole system
`j_ver`                     |  kg/m^2 s | Mean value of data in `jv_vector.txt`
`j_hor`                     |  kg/m^2 s | Mean value of data in `jh_vector.txt`
`Time`                      |   N/A     | Time to run the software [Days/Hours/Minutes/Seconds]
