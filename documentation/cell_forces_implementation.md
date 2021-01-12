# Implementing custom cell-generated forces

We aim to implement a function that, at each time step, generates a new force value based on the velocity distributions found in Plou et al. (2018).

We will be using inverse transform sampling to generate new values from a uniformly distributed interval. At each velocity update (this function will be integrated in the velocity update function), we will be computing a random value between 0 and 1 and transforming this value according to a previously defined polynomial function.

:pencil: **Code**

```C++
double locomotive_forces_generator( )
{
    // random number generator to define cell locomotive forces
    // based on anempirically obtained velocity distribution
    double random_value, force_value;
    random_value = UniformRandom();
    force_value =
    0.1157 * pow(random_value, 3) + 0.2422 * pow(random_value, 2) +
    0.0053 * random_value + 0.0048;
    force_value = 13.5*force_value;
    return force_value;
}

```

## References

Plou, J., et al. "From individual to collective 3D cancer dissemination: roles of collagen concentration and TGF-β." Scientific reports 8.1 (2018): 1-14
