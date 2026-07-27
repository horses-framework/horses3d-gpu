title: Postprocessing 
---

[TOC]

For postprocessing the Simulation Results

## Visualization with Tecplot Format: *horses2plt*


`HORSES3D` provides a script for converting the native binary solution files (*.hsol) into tecplot ASCII format (*.tec), which can be visualized in Pareview or Tecplot. It can also export the solution to the more recent VTKHDF format; however, note that this feature does **not** export boundary information or mesh files. Usage:

```bash
	$ horses2plt SolutionFile.hsol MeshFile.hmesh <<Options>>
```

The options comprise following flags:

| Flag                | Description                                                                                                      | Default value |
|---------------------|------------------------------------------------------------------------------------------------------------------|---------------|
| --output-order=    | *INTEGER*: Output order nodes. The solution is interpolated into the desired number of points.                  | Not Present   |
| --output-basis=    | *CHARACTER*: Either *Homogeneous* (for equispaced nodes, or *Gauss*.)                                            | *Gauss*       |
| --output-mode=     | *CHARACTER*: Either *multizone* or *FE*. The option *multizone* generates a Tecplot zone for each element. The option *FE* generates only one Tecplot zone for the fluid and one for each boundary (if *--boundary-file* is defined). Each subcell is mapped as a linear finite element. This format is faster to read by Paraview and Tecplot. | *multizone*   |
| --output-variables=| *CHARACTER*: Output variables separated by commas. A complete description can be found in Section 2.             | Q             |
| --dimensionless    | Specifies that the output quantities must be dimensionless.                                                     | Not Present   |
| --partition-file=  | *CHARACTER*: Specifies the path to the partition file (*.pmesh) to export the MPI ranks of the simulation.     | Not Present   |
| --boundary-file=   | *CHARACTER*: Specifies the path to the boundary mesh file (*.bmesh) to export the surfaces as additional zones of the Tecplot file. | Not Present   |
| --output-type=     | *CHARACTER*: Specifies the type of output file: *tecplot* or *vtkhdf*.                                          | *tecplot*     |

* *Homogeneous* when *--output-order* is specified

Additionally, depending on the type of solution file, the user can specify additional options.

## Solution Files (*.hsol) 

For standard solution files, the user can specify which variables they want to be exported to the Tecplot file with the flag **-{}-output-variables=*.
The options are:

<div class="multicols" style="column-count: 5;">
  <ul>
    <li>\(Q\) (default)</li>
    <li>\(\rho\)</li>
    <li>\(u\)</li>
    <li>\(v\)</li>
    <li>\(w\)</li>
    <li>\(p\)</li>
    <li>\(T\)</li>
    <li>\(Mach\)</li>
    <li>\(s\)</li>
    <li>\(Vabs\)</li>
    <li>\(V\)</li>
    <li>\(Ht\)</li>
    <li>\(rhou\)</li>
    <li>\(rhov\)</li>
    <li>\(rhow\)</li>
    <li>\(rhoe\)</li>
    <li>\(c\)</li>
    <li>\(Nxi\)</li>
    <li>\(Neta\)</li>
    <li>\(Nzeta\)</li>
    <li>\(Nav\)</li>
    <li>\(N\)</li>
    <li>\(Ax\_Xi\)</li>
    <li>\(Ax\_Eta\)</li>
    <li>\(Ax\_Zeta\)</li>
    <li>\(ThreeAxes\)</li>
    <li>\(Axes\)</li>
    <li>\(mpi\_rank\)</li>
    <li>\(eID\)</li>
    <li>\(gradV\)</li>
    <li>\(u\_x\)</li>
    <li>\(v\_x\)</li>
    <li>\(w\_x\)</li>
    <li>\(u\_y\)</li>
    <li>\(v\_y\)</li>
    <li>\(w\_y\)</li>
    <li>\(u\_z\)</li>
    <li>\(v\_z\)</li>
    <li>\(w\_z\)</li>
    <li>\(c\_x\)</li>
    <li>\(c\_y\)</li>
    <li>\(c\_z\)</li>
    <li>\(\omega\)</li>
    <li>\(\omega\_x\)</li>
    <li>\(\omega\_y\)</li>
    <li>\(\omega\_z\)</li>
    <li>\(\omega\_abs\)</li>
    <li>\(Qcrit\)</li>
  </ul>
</div>


## Surface Solution Files (*.surf.hsol)

`HORSES3D` can save the solution restricted to a set of boundaries (instead of the full volume) at a fixed time interval, while the simulation is running. These *surface solution files* are the recommended way of postprocessing wall quantities such as the friction velocity, the wall shear stress, the skin friction coefficient \(C_f\) or \(y^+\), since they are much smaller than full volume snapshots and can be saved much more frequently.

### Generating surface solution files

Surface saving is configured entirely from the **simulation** control file (the one used to run the solver, not the one used by *horses2plt*). The relevant keywords are:

| Keyword                        | Description                                                                                                                     | Default value |
|---------------------------------|----------------------------------------------------------------------------------------------------------------------------------|---------------|
| `boundaries to save`            | *ARRAY*: List of boundary marker names for which a surface solution file is generated, e.g. `[cylinder]`.                        | Not present   |
| `surface save timestep`         | *REAL*: Simulation-time interval between surface autosaves.                                                                      | Not present   |
| `surface save gradients`        | *LOGICAL*: Additionally save the velocity gradients (\(\partial u_i/\partial x_j\)) on the surface.                              | `.false.`     |
| `surface save utau`             | *LOGICAL*: Save the (scalar) friction velocity \(u_\tau\) on no-slip walls.                                                      | `.false.`     |
| `surface save utau vector`      | *LOGICAL*: Save the friction-velocity **vector** (Cartesian components \(u_{\tau,x}, u_{\tau,y}, u_{\tau,z}\)) on no-slip walls. | `.false.`     |
| `surface save turbulent`        | *LOGICAL*: Save the wall viscosity (`mu_ns`) and the wall-normal distance (`wall_distance`), needed for \(y^+\).                 | `.false.`     |

Only no-slip wall boundaries actually receive the `utau`/`utau vector`/`turbulent` data; the other flags apply to any saved boundary (walls or slices). `surface save utau` and `surface save utau vector` are independent and can both be active at the same time. Example:

```
boundaries to save       = [cylinder]
surface save timestep    = 1e-3
surface save gradients   = .false.
surface save turbulent   = .true.
surface save utau        = .true.
surface save utau vector = .true.
```

**Important:** these flags only control what gets *written* to the `.surf.hsol` file at solve time. The corresponding `additional variables` flags used by *horses2plt* (see below) only control what it *expects to read* — they do not create data that was not saved. Requesting a variable in *horses2plt* that was not actually saved by the solver will misalign the reads for every following face and typically causes a crash or an explicit "Array found in file dimensions does not match..." error. Always keep the solver's `surface save ...` flags and *horses2plt*'s `additional variables` in sync, and rebuild/rerun the **solver** binary (not just *horses2plt*) after any code change that affects these fields, since the writing code lives in the shared mesh library.

### Postprocessing surface solution files

Surface files are postprocessed with the same *horses2plt* tool, pointing it at the surface mesh/solution instead of the volume ones:

```
hmesh file = "MESH/Cylinder_bc001_0000000000.surf.hmesh"
hsol file  = "RESULTS/surfaces/Cylinder_bc001_0000000906.surf.hsol"
```

Wall/turbulence-related quantities are not read by default: they must be requested with the `additional variables` keyword, matching what the solver actually saved:

| `additional variables` token | Enables reading                                   | Requires the solver flag(s)                          |
|-------------------------------|----------------------------------------------------|-------------------------------------------------------|
| `u_tau`                       | Scalar friction velocity                           | `surface save utau = .true.`                           |
| `u_tau_vector`                 | Friction-velocity vector (Cartesian components)    | `surface save utau vector = .true.`                    |
| `turb`                         | Wall viscosity + wall distance (needed for \(y^+\)) | `surface save turbulent = .true.`                      |
| `les`                          | Sub-grid viscosity (`mu_sgs`)                       | (LES/SGS output, volume files)                          |

```
additional variables = "[u_tau,u_tau_vector,turb]"
output variables = Cp, Cf, yplus, wall_distance, u_tau
```

### Wall/turbulence-related output variables

<div class="multicols" style="column-count: 2;">
  <ul>
    <li>\(u\_tau\)</li>
    <li>\(u\_tau\_x\)</li>
    <li>\(u\_tau\_y\)</li>
    <li>\(u\_tau\_z\)</li>
    <li>\(wall\_shear\)</li>
    <li>\(Cf\)</li>
    <li>\(yplus\)</li>
    <li>\(wall\_distance\)</li>
    <li>\(mu\_ns\)</li>
    <li>\(mu\_sgs\)</li>
  </ul>
</div>

- **`u_tau`**: friction velocity magnitude, \(u_\tau = \sqrt{|\tau_w|/\rho}\) (velocity units). Requires the `u_tau` additional variable.
- **`u_tau_x`, `u_tau_y`, `u_tau_z`**: Cartesian components of the friction-velocity vector, same velocity units and sign convention as the local wall shear stress direction (unlike the scalar `u_tau`, which is always reported as a magnitude). By construction, \(\sqrt{u\_tau\_x^2+u\_tau\_y^2+u\_tau\_z^2} = u\_tau\). Requires the `u_tau_vector` additional variable.
- **Requesting `u_tau` as an output variable is context-dependent**: asking for `output variables = ... u_tau ...` does **not** always produce a single column — it expands automatically depending on which additional variables were activated:
    - only `u_tau` active &rarr; 1 column: `u_tau`.
    - only `u_tau_vector` active &rarr; 3 columns: `u_tau_x, u_tau_y, u_tau_z`.
    - both active &rarr; 4 columns: `u_tau, u_tau_x, u_tau_y, u_tau_z`.

  There is no separate `u_tau_vector` output variable to request; the single `u_tau` key covers every combination.
- **`wall_shear`**: wall shear stress \(\tau_w = \rho\,u_\tau^2\,\text{sign}(u_\tau)\) (stress units). Requires `u_tau`.
- **`Cf`**: skin friction coefficient, \(C_f = 2\,\tau_w\) in the code's reference-normalized nondimensional units, equivalent to the standard definition \(C_f = \tau_w/(\tfrac{1}{2}\rho_\infty U_\infty^2)\). Requires `u_tau`. This is **not** algebraically the same as reconstructing `Cf` from `u_tau`/`u_tau_x/y/z` divided by a reference velocity (that reintroduces the *local* density instead of the reference one — see the discussion below).
- **`yplus`**: wall-normal distance in wall units. Requires `u_tau` and `turb`.
- **`wall_distance`**: distance from the first fluid solution point to the wall. Requires `turb`.
- **`mu_ns`**: molecular (+ SGS, if LES) viscosity at the wall. Requires `turb` or `les`.

#### `Cf` vs. a manual `u_tau`-based reconstruction

Because \(u_\tau\) is defined using the **local** density, \(u_\tau=\sqrt{|\tau_w|/\rho_{local}}\), reconstructing a skin friction coefficient by hand as `2*(u_tau/U_ref)^2` (e.g. in ParaView, from `u_tau` or `sqrt(u_tau_x^2+u_tau_y^2+u_tau_z^2)`) computes \(2\,\tau_w/\rho_{local}\), **not** the standard \(C_f = 2\,\tau_w/\rho_\infty\) that the `Cf` output variable already provides directly. The two only coincide where \(\rho_{local}\approx\rho_\infty\), and can differ substantially in regions with strong local compressibility effects (e.g. across a shock, or a stagnation point). Prefer the `Cf` output variable over a manual reconstruction unless you specifically need it normalized by the local density (e.g. to compare against another code that itself uses that convention).

## Statistics Files (*.stats.hsol)
Statistics files can generate the standard variables as well as the following variables (being \(S_{ij}\) the components of the Reynolds Stress tensor):

<div class="multicols" style="column-count: 3;">
  <ul>
    <li>\(umean\)</li>
    <li>\(vmean\)</li>
    <li>\(wmean\)</li>
    <li>\(S_{xx}\)</li>
    <li>\(S_{yy}\)</li>
    <li>\(S_{zz}\)</li>
    <li>\(S_{xy}\)</li>
    <li>\(S_{xz}\)</li>
    <li>\(S_{yz}\)</li>
  </ul>
</div>



## Extract geometry
Under construction.

## Merge statistics tool

Tool to merge several statistics files. The usage is the following:

```bash
	$ horses.mergeStats *.hsol --initial-iteration=INTEGER --file-name=CHARACTER
```

Some remarks:

- Only usable with statistics files that are obtained with the "reset interval" keyword and/or with individual consecutive simulations.
- Only constant time-stepping is supported.
- If the hsol files have the gradients, the following flag must be used
```bash
	$ --has-gradients
```
- Dynamic p-adaptation is currently not supported.

## Mapping result to different mesh 
<a name="MaptomeshKey"></a>
`HORSES3D` addons, *horsesConverter*, has a capability to map result into different mesh file, with both have a consistent geometry. This is done by performing interpolation with the polynomial inside each element for each node point of the new mesh. The type of node quadrature will follow the quadrature defined in the .hmesh file with selected polynomial order in the control file. A control input file is required and must has name *horsesConverter.convert*. The template of control input file will be generated by default when executing *./horsesConverter* in a directory without the control file. Error message is given when at least one node point of the new mesh is not within any element of the old mesh. After completion, a new result file is generated and named *Result\_interpolation.hsol*. The required keywords in the control file are described in the table below. Command to execute:
```bash
	$ ./horsesConverter
```

| Keyword              | Description                                                        | Default value |
|----------------------|--------------------------------------------------------------------|---------------|
| Task                 | *meshInterpolation*                                               |               |
| Mesh Filename 1      | Location of the origin mesh (*.hmesh)                             |               |
| Boundary Filename 1  | Location of the origin boundary mesh (*.bmesh)                    |               |
| Result 1             | Location of the solution file with origin mesh (*.hsol)           |               |
| Mesh Filename 2      | Location of the target mesh (*.hmesh)                             |               |
| Boundary Filename 2  | Location of the target boundary mesh (*.bmesh)                    |               |
| Polynomial Order     | Polynomial order of the target mesh                                | (1, 1, 1)     |


## Generate OpenFOAM mesh 
<a name="GenerateOpenFOAMmeshKey"></a>
Another functionality of *horsesConverter* addons is to convert the mesh files, (\*.hmesh) and (\*.bmesh), into OpenFOAM format, the *polyMesh* folder. Each element is discretised into \(n_x \times n_y \times n_z\) cells distributed as Gauss-Lobatto nodes. The number of division of each element, (\(n_x\), \(n_y\), and \(n_z\)), is required in the control file, see [previous section](GenerateOpenFOAMmeshKey). After completion, a folder named `foamFiles` is generated. OpenFOAM mesh files, i.e. *points*, *faces*, *owner*, *neighbour*, and *boundary*, are located within `foamFiles/constant/polyMesh`. The required keywords in the control file are described in the table below Command to execute:

```bash
	$ ./horsesConverter
```

| Keyword              | Description                                                        | Default value |
|----------------------|--------------------------------------------------------------------|---------------|
| Task                 | *horsesMesh2OF*                                                   |               |
| Mesh Filename 1      | Location of the origin mesh (*.hmesh)                             |               |
| Boundary Filename 1  | Location of the origin boundary mesh (*.bmesh)                    |               |
| Polynomial Order     | Number of division of each element (\(n_x\), \(n_y\), and \(n_z\))      | (1, 1, 1)     |

NOTE: Before running the mesh in the OpenFOAM environment, the type of boundaries inside the boundary file needs to be adjusted according to the actual type (*patch*, *wall*, and *symmetry*).


## Generate HORSES3D solution file from OpenFOAM result
`HORSES3D` provides a capability to convert OpenFOAM result into `HORSES3D` solution file (\*.hsol). The mesh of the OpenFOAM result must be generated by converting HORSES3D mesh files, see the [previous section](GenerateOpenFOAMmeshKey). Beforehand, the OpenFOAM result must be converted into VTK format(*.vtk). This not only allows the result to be in the single file but also converts cell data into point data. In the OpenFOAM environment, the command for this conversion:  
```cpp
	$ foamToVTK -fields "(U p T rho)" -ascii -latestTime
```
The necessary file (.vtk) required for the control file input is inside VTK folder, see [previous section](GenerateOpenFOAMmeshKey) for the control file template. The `HORSES3D` solution file is named `Result\_OF.hsol`. The required keywords in the control file are described in the table below. Command to execute:

```bash
	$ ./horsesConverter
```

| Keyword                         | Description                                               | Default value |
|---------------------------------|-----------------------------------------------------------|---------------|
| Task                            | *OF2Horses*                                               |               |
| Mesh Filename 1                 | Location of the origin mesh (*.hmesh)                    |               |
| Boundary Filename 1             | Location of the origin boundary mesh (*.bmesh)           |               |
| Polynomial Order                | Polynomial order of the solution file (.hsol)             | (1, 1, 1)     |
| VTK file                        | Location of VTK file (.vtk)                               |               |
| Reynolds Number                 | Reynolds Number/m of the solution -- \(L_{ref}\)=1.0m       |               |
| Mach Number                     | Mach Number of the solution                               |               |
| Reference pressure (Pa)         | Reference Pressure                                        | 101325        |
| Reference temperature (K)       | Reference Temperature                                     | 288.889       |

