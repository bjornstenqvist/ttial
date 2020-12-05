# ...The Tea Is A Lie
### ...The Tea Is False
#### ...The Tea = 0
##### ...d Tea = 0
###### ...d t = 0
...steady state

This is a C++ software for calculating membrane concentration-profiles, fluxes, permeabilities, and so forth in steady-state based on Mesh Analysis. The software is built around a 2D rectangular membrane which is laterally periodic and therefore infinite, suitable for for example brick and mortar systems. The program is fundamentally constructed to model mass diffusion however it can also be used to model temperature flux and change diffusion. The major restriction of the scheme is that it assumes one constant concentration profile at the upper boundary of the system while utilizing a zink, i.e. zero concnentration, at the lower boundary. For thermal/electric diffusion the temperature/electric potential is constant at those boundaries.

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

### Input/Output data

In the following we define the diffusion coefficient as D=U*R*T where U is the mobility, R the gas constant, and T the temperature. The solubility is defined by S = c0 / gamma / exp(mu0/RT) where c0 is the standrad concentration, gamma the activity coefficient, and mu0 the standard chemical potential. For charge diffusion `c_out` divited by `S_out` is the electric potential difference over the system and the resistances are the products between the solubilities `S` and diffusion coefficients `D`. For heat diffusion `c_out` divited by `S_out` is the temperature difference over the system and the thermal conductivities are the products between the solubilities `S` and diffusion coefficients `D`.

Input parameter    | Description [unit]
------------------ | -------------------
`model_nbr`	   | model number, 0 gives Brick and Mortar [integer]
`d`		   | width of bricks [double]
`s`		   | horizontal spacing between bricks [double]
`g`		   | vertical spacing between bricks [double]
`t`		   | brick thickness [double]
`N`		   | nbr of layers of bricks [integer]
`omega`		   | offset ratio, negative gives random [double]
`c_out`		   | concentration outside system [double]
`S_out`		   | solubility outside system [double]
`S_mv`		   | mortar solubility (verticle) [double]
`S_mh`		   | mortar solubility (horizontal) [double]
`S_bv`		   | brick solubility (verticle) [double]
`S_bh`		   | brick solubility (horizontal) [double]
`D_mv`		   | mortar diffusion coefficient (verticle) [double]
`D_mh`		   | mortar diffusion coefficient (horizontal) [double]
`D_bv`		   | brick diffusion coefficient (verticle) [double]
`D_bh`		   | brick diffusion coefficient (horizontal) [double]
`seed`		   | seed, negative gives random [integer]
`height`	   | height of system, only relevant when loading an external mesh [double]
`width`		   | width of system, only relevant when loading an external mesh [double]
`load_external`	   | load external mesh, true/false [bool]
`C`		   | number of columns of nodes [integer]
`R`		   | number of rows of nodes [integer]

Output files                | Description
--------------------------- | -------------
`conc_hor_matrix.txt`       | Concentration in the horizontal cells
`conc_ver_matrix.txt`       | Concentration in the vertical cells
`sh_matrix.txt`             | Solubility in the horizontal cells
`sv_matrix.txt`             | Solubility in the vertical cells
`j_hor_matrix.txt`          | Flux in the horizontal cells
`j_ver_matrix.txt`          | Flux in the vertical cells
`j_hor_vector.txt`          | Mean of amplitude of flux in each row of the horizontal cells
`j_ver_vector.txt`          | Sum of flux in each row of the vertical cells
`Dh_matrix.txt`             | Diffusion coefficients in the horizontal cells
`Dv_matrix.txt`             | Diffusion coefficients in the vertical cells
`V_matrix.txt`              | Potential in nodes in between cells
`output.txt`                | Miscellaneous data

Data in `output.txt`        | Description
--------------------------- | -------------
`model`                     | Name of used model
`K_{M_ver/out}`             | Partition coefficient between mortar (vertical) and outside
`K_{B_ver/out}`             | Partition coefficient between bricks (vertical) and outside
`K_{M_ver/B_ver}`           | Partition coefficient between mortar and bricks (both vertical) 
`K_{M_hor/B_hor}`           | Partition coefficient between mortar and bricks (both horizontal) 
`I_bot_sum`                 | Flux out of the system
`I_top_sum`                 | Flux into the system
`Reff`                      | Effective resistance of the whole system
`j_ver`                     | Mean value of data in `jv_vector.txt`
`j_hor`                     | Mean value of data in `jh_vector.txt`
`Time`                      | Time to run the software [Days/Hours/Minutes/Seconds]
