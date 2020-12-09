
model BIOMD178 "Lebeda2008 - BoTN Paralysis (4 step model)"

    parameter Real kT = 0.141;
    parameter Real kL = 0.013;
    parameter Real kB = 0.058;
    parameter Real kS = 0.00015;

    Real tension(start=0.0);

    Real BoNT;
    Real bulk;
    Real free;
    Real bound;
    Real translocate;
    Real lytic;

initial equation
    BoNT = 0.0;
    bulk = 1.0;
    free = 0.0;
    bound = 0.0;
    translocate = 0.0;
    lytic = 0.0;

equation
    tension = 1 - lytic
    BoNT = bulk + free
    der(bulk) =  - (kS * bulk)
    der(free) = (kS * bulk) - (kB * free)
    der(bound) = (kB * free) - (kT * bound)
    der(translocate) = (kT * bound) - (kL * translocate)
    der(lytic) = (kL * translocate) 

end BIOMD178;
