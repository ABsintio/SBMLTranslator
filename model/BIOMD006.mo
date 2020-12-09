
model m "Tyson1991 - Cell Cycle 2 var"

    parameter Real kappa = 0.015;
    parameter Real k6 = 1.0;
    parameter Real k4 = 180.0;
    parameter Real k4prime = 0.018;

    Real alpha(start=0.0);

    Real EmptySet;
    Real u;
    Real z;
    Real v;

initial equation
    EmptySet = 0.0;
    u = 0.0;
    z = 0.0;
    v = 0.0;

equation
    z = v - u
    alpha = k4prime / k4
    der(u) = k4 * (v - u) * (alpha + pow(u, 2)) - k6 * u
    der(v) = kappa - k6 * u
    der(EmptySet) = 0.0

end m;
