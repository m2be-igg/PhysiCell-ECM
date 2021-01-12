# Introducing the ECM and its effect on cell motility

## The ECM as a non-diffusing substance

### Defining the ECM

The ECM is introduced in the `PhysiCell_settings.xml` (see the ./config/ folder) file as a non-diffusing substance. We define the ECM by its collagen density. Thus, in the XML file we will introduce a new substance with a defined concentration. Later, we will use this value to compute the mechanical properties of the matrix.

:pencil: **Code**
Defining the ECM as a matrix with a collagen density of 2.5 mg/mL.

```xml
<variable name="ECM" units="mg/mL" ID="1">
    <physical_parameter_set>
        <diffusion_coefficient units="micron^2/min">
            0.0
        </diffusion_coefficient>
        <decay_rate units="1/min">
            0.0
        </decay_rate>
    </physical_parameter_set>
    <initial_condition units="mg/mL">
        2.5
    </initial_condition>
    <Dirichlet_boundary_condition units="mg/mL" enabled="false">
        2.5
    </Dirichlet_boundary_condition>
</variable>
```

:bulb: **Future implementations**

- Non-uniform initial distribution
- Decay rate

## The ECM as a regulator of cell motility

We introduce the effect of the ECM on cell motility in the `custom.cpp` file. Our goal is to take into account the dynamic viscosity of the matrix when cell velocity is computed. To do so, we need to adapt the original cell velocity update function (`standard_update_cell_velocity`) and divide its output by the dynamic viscosity of the matrix.

To compute the dynamic viscosity value, we sample the ECM concentration at each voxel. Based on this concentration value, we define the dynamic viscosity value, taking into consideration the experimental measurements made in Valero et al. (2018).

:pencil: **Code**

```c++
void drag_update_cell_velocity( Cell* pCell, Phenotype& phenotype, double dt ) {
    // sample ECM
    int ECM_density_index = microenvironment.find_density_index( "ECM" );
    double ECM_density = pCell->nearest_density_vector()[ECM_density_index];
    double dyn_viscosity;
    // get viscosity based on concentration
    if(ECM_density == 2.5) {
        dyn_viscosity = 7.96;
    }
    else if(ECM_density == 4.0){
        dyn_viscosity = 18.42;
    }
    else if(ECM_density == 6.0) {
        dyn_viscosity = 39.15;
    }
    // update the speed value
    pCell->phenotype.motility.migration_speed = 
        locomotive_forces_generator();
    // update velocity
    standard_update_cell_velocity(pCell, phenotype, dt);
    // include the 1/vu (1/ECM density) term to consider friction
    pCell->velocity /= dyn_viscosity;
    return;
}

```

The initial part of this code can be disregarded for a uniform distribution of collagen, in which there is no decay nor decay. Instead, the dynamic viscosity could be defined at the beginning of the simulaitons and be kept constant through time. However, we have written this function to be able to extended for non-uniformly distributed domains.

:bulb: **Future implementations**

- Continuum function to relate collagen density and dynamic viscosity values

## References

Valero, Clara, et al. "Combined experimental and computational characterization of crosslinked collagen-based hydrogels." PLoS One 13.4 (2018): e0195820.
