[![Testsuite](https://github.com/ManelSoria/HGSpy/actions/workflows/run_examples.yml/badge.svg)](https://github.com/CalebFuster/HGSpy/actions)
[![License](https://img.shields.io/badge/license-GPL3-orange)](https://opensource.org/license/gpl-3-0/)

# HGSpy

This package aims to provide an **open source** and **comprehensive** approach to combustion
problems and further isentropic expansion using the *<span style="color:Aquamarine">NASA polynomials*
To do it, this module is separated into 3 parts: the [thermodynamic properties](#prop),
the [Algorithms](#alg) and [Other functionalities](#other).

This is the python implementation of the original HGS code in MATLAB that can be found [in Github](https://github.com/ManelSoria/HGS) or [MATLAB File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/40862-hgs-chemical-equation-solver)

The original script is developed/maintained by 
[C. Fuster](https://www.linkedin.com/in/caleb-fuster-b6a566182/), 
[M. Soria](https://directori.upc.edu/directori/dadesPersona.jsp?id=1003031), 
[A. Miró](https://directori.upc.edu/directori/dadesPersona.jsp?id=1126813), 
et al.

<span style="color:Aquamarine">NASA polynomials</span> are downloaded from [Burcat mirrors link](http://garfield.chem.elte.hu/Burcat/THERM.DAT) check his [documentation](http://garfield.chem.elte.hu/Burcat/Archives/2005.pdf).

There is no more documentation in this package that the comments on the functions, but you can check the master thesis developed with the Matlab code for more info.
Request it to author of this package or download it from [UPCcommons](https://upcommons.upc.edu/discover?filtertype_1=author&filter_relational_operator_1=contains&filter_1=Caleb+Fuster&submit_apply_filter=)
when the institution publishes it.

## Deployment

A _Makefile_ is provided within the tool to automate the installation for easiness of use for the user. To install the tool simply create a virtual environment as stated below or use the system Python. Once this is done simply type:
```bash
make
```
This will install all the requirements and install the package to your active python. To uninstall simply use
```bash
make uninstall
```

The previous operations can be done one step at a time using
```bash
make requirements
```
to install all the requirements;
```bash
make python
```
to compile and;
```bash
make install
```
to install the tool.

### Virtual environment

The package can be installed in a Python virtual environement to avoid messing with the system Python installation.
Next, we will use [Conda](https://docs.conda.io/projects/conda/en/latest/index.html) for this purpose.
Assuming that Conda is already installed, we can create a virtual environment with a specific python version and name (`my_env`) using
```bash
conda create -n my_env python=3.8
```
The environment is placed in `~/.conda/envs/my_env`.
Next we activate it be able to install packages using `conda` itself or another Python package manager in the environment directory:
```bash
conda activate my_env
```
Then just follow the instructions as stated above.


## Modules

---
<h2 id="prop"> Thermodynamic properties </h2>

Thermodynamic properties are accessible  using 
<span style="color:tomato">*HGS.single()*</span> or <span style="color:tomato">*HGS.prop()*</span>
(see their own documentation using **<span style="color:grey">help()</span>**)


<details>
<summary><span style="font-weight: bold;color:grey">help</span>(<span style="color:tomato">fun</span>)</summary>
<details id="hgssingle">
    <summary><span style="color:tomato">HGS.single</span></summary>

    *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

    res = hgs_single(species, prop, T, P)

    *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

    hgs_single returns the property of a species

    *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
    Inputs:
    -----------------------------------------------------------------------------
    species --> String or numbers of species
    prop --> Property requested (see below)
    T --> [K] Temperature
    P --> [bar] Pressure

    *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
    Outputs:
    -----------------------------------------------------------------------------
    res --> Property result
          mm [g/mol]
          cp [kJ/(mol*K)]
          cv [kJ/(mol*K)]
          h [kJ/mol]
          s [kJ/(mol*K)]
          g [kJ/mol]

    *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
    * Matlab HGS 2.0
    * By Caleb Fuster, Manel Soria and Arnau Miró
    * ESEIAAT UPC
</details>
<details id="hgsprop">
    <summary><span style="color:tomato">HGS.prop</span></summary>

    *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

    var = hgs_prop(species, n, T, P, *args)

    *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

    hgs_prop returns the properties of the mixture of gasses

    *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
    Inputs:
    -----------------------------------------------------------------------------
    species --> String or numbers of species
    prop --> Property requested (see below)
    T --> [K] Temperature
    P --> [bar] Pressure
    *args --> Expected return: 'Mm' 'Cp' 'Cv' 'H' 'S' 'G' 'Rg' 'gamma' 'a'
                               If it is empty, all the properties will be
                               return

    *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
    Outputs:
    -----------------------------------------------------------------------------
    var --> Property result
          mm [g/mol]
          cp [kJ/K]
          cv [kJ/K]
          h [kJ]
          s [kJ/K]
          g [kJ]
          Rg [kJ/(kg*K)]
          gamma
          a [m/s]

    *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
    * Python HGS 1.0 from Matlab HGS 2.0
    * By Caleb Fuster, Manel Soria and Arnau Miró
    * ESEIAAT UPC
</details>
</details>



---
<h2 id="alg"> Algorithms </h2>

This module includes an equilibrium algorithm (<span style="color:tomato">*HGS.eq()*</span>), 
combustion algorithm (<span style="color:tomato">*HGS.tp()*</span>) and 
isentropic expansion algorithm (<span style="color:tomato">*HGS.isentropic()*</span>)
(see their own documentation using <span style="color:grey">**help()**</span>)


<details>
<summary><span style="font-weight: bold;color:grey">help</span>(<span style="color:tomato">fun</span>)</summary>
<details id="hgseq">
    <summary><span style="color:tomato">HGS.eq</span></summary>

     *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

    species, n, Gmin = hgs_eq(species, n0, T, P, **kwargs)

    *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

    hgs_eq calculates the species mols equilibrium at a certain temperature
    and pressure

    *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
    Inputs:
    -----------------------------------------------------------------------------
    species --> String or numbers of species
    n0 --> [mol] Initial mixture
    T --> [K] Temperature. Could be a single value or an array.
    P --> [bar] Pressure
    **kwargs --> opti_eq= Options for the minimize Scipy function

    *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
    Outputs:
    -----------------------------------------------------------------------------
    species --> Species
    n --> [mol] Final mixture
    Gmin --> [kJ] Minimum Gibbs free energy

    *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
    * Python HGS 1.0 from Matlab HGS 2.0
    * By Caleb Fuster, Manel Soria and Arnau Miró
    * ESEIAAT UPC
</details>
<details id="hgstp">
    <summary><span style="color:tomato">HGS.tp</span></summary>

    *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

    Tp, n, species, flag = hgs_tp(species, n0, tipo, V0, P, **kwargs)

    *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

    hgs_tp calculates the reaction temperature considering dissociation,
    and the products composition in equilibrium

    *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
    Inputs:
    -----------------------------------------------------------------------------
    species --> String or numbers of species
    n0 --> [mols] Number of mols of each species
    tipo --> Entry type that defines the state of the input.
             It can be 'T' or 'H'
    V0 --> Entry that should be for type:'T'   V0=T [K] input temperature
                                         'H'   V0=H [kJ] input enthalpy
    P --> [bar] Mixture pressure
    **kwargs --> opti_eq= Options for the minimize Scipy function
                 opt_sec= Dictionary with the options for the secant method.
                        "xmin": [K] Temperature minimum for the solver;
                        "xmax" [K] Temperature maximum for the solver;
                        "maxiter" Max iterations for the solver;
                        "epsx" Diferential T where the solver reachs the solution;
                        "epsy" Diferential S where the solver reachs the solution;
                        "fchange" T difference where secant method is
                                 changed by bisection method;
                        "tipo" Select between: 'Frozen' for frozen flow
                                              'Shifting' for shifting flow
                        "info" Detailed info == 1; No info == 0.
                        "dTp" Improve the velocity with the approximation of
                              parabola. +- dTp
                        opt_sec = {"xmin": 300, "xmax": 4000, "maxiter": 200,
                                   "epsx": 0.1, "epsy": 1, "tipo": "Shifting",
                                   "fchange": 5, "info": 0, "dTp": 100}

    *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
    Outputs:
    -----------------------------------------------------------------------------
    Tp --> [K] Final temperature
    n --> [mol] Final mixture
    species --> String or numbers of species
    flag --> Solver error detection:
                  1  Solver has reached the solution
                 -1  Solver failed. Maximum iterations
                 -2  Solver failed. Initial sign change not found

    *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
    * Python HGS 1.0 from Matlab HGS 2.0
    * By Caleb Fuster, Manel Soria and Arnau Miró
    * ESEIAAT UPC
</details>
<details id="hgsisentropic">
    <summary><span style="color:tomato">HGS.isentropic</span></summary>

    *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

    Tp, n, species, v2, M2, flag = hgs_isentropic(species, n0, T0, P0, P1, **kwargs)

    *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

    hgs_tp calculates the outlet variables for an isentropic expansion

    *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
    Inputs:
    -------------------------------------------------------------------------------
    species --> String or numbers of species
    n0 --> [mols] Number of mols of each species
    T0 --> [K] Initial temperature
    P0 --> [bar] Inlet pressure
    P1 --> [bar] Exit pressure
    **kwargs --> opti_eq= Options for the minimize Scipy function
                 opt_sec= Dictionary with the options for the secant method.
                        "xmin": [K] Temperature minimum for the solver;
                        "xmax" [K] Temperature maximum for the solver;
                        "maxiter" Max iterations for the solver;
                        "epsx" Diferential T where the solver reachs the solution;
                        "epsy" Diferential S where the solver reachs the solution;
                        "fchange" T difference where secant method is
                                 changed by bisection method;
                        "tipo" Select between: 'Frozen' for frozen flow
                                              'Shifting' for shifting flow
                        "info" Detailed info == 1; No info == 0.
                        "dTp" Improve the velocity with the approximation of
                              parabola. +- dTp
                        opt_sec = {"xmin": 300, "xmax": 4000, "maxiter": 200,
                                   "epsx": 0.1, "epsy": 1, "tipo": "Shifting",
                                   "fchange": 5, "info": 0, "dTp": 100}

    *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
    Outputs:
    -----------------------------------------------------------------------------
    Tp --> [K] Exit temperature
    n --> [mols] Species resultant mols
    species --> String or numbers of species
    v2 --> [m/s] Velocity of the mixture
    M2 --> [Mach] Mach of the mixture
    flag --> Solver error detection:
                  1  Solver has reached the solution
                 -1  Solver failed. Maximum iterations
                 -2  Solver failed. Initial sign change not found

    *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
    * Python HGS 1.0 from Matlab HGS 2.0
    * By Caleb Fuster, Manel Soria and Arnau Miró
    * ESEIAAT UPC
</details>
</details>

---
<h2 id="other"> Other functionalities</h2>

This module includes other functionalities to facilitate the user life. 
You can:

Update the database:

```python
import HGSpy as HGS
DOWNLOAD = False # Whether to rebuild or just download the database
INFO     = True  # Print information when downloading
db = HGS.data.download(download=DOWNLOAD,info=INFO)
db.save(HGS.HGSDATA)
```

Create/Erase your own mixtures:

```python
import HGSpy as HGS
db = HGS.HGSData.load()
db.add("NAME", ["species1", "species2",...], [20, 80,...])
db.subt("NAME")
db.save(HGS.HGSDATA)
```

And search and print info from a species:

```python
import HGSpy as HGS
HGS.find("NAME")
HGS.print_info("NAME")
```

