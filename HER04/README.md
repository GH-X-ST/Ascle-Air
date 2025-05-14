Systems, Propulsion and Mission Equipment

# Breakdown of Files:

## conceptual_weights.mlx

MATLAB Live script to combine and calculate total weight, rotorcraft CG and moments of inertia.<br>

## mass_input.xlsx

Excel document for component or system weights to be input as point masses with displacements from the reference point of the rotorcraft nose.<br>

### Required Inputs:

m: mass in kilograms<br>

x: x-displacement from reference point in metres<br>

y: y-displacement from reference point in metres<br>

z: z-displacement from reference point in metres<br>

i_xx, i_yy, i_zz: Moments of Inertia in the x, y and z axes about the component's centre of mass, if the component is to be modelled as a point mass, leave zero<br>

i_xy, i_xz, i_zy: Product Moments of Inertia about the component's centre of mass, if the component is to be modelled as a point mass, leave zero<br>