
model BIOMD040 ""

    parameter Real f = 1.0;
    parameter Real k1 = 1.34;
    parameter Real k2 = 1600000000.0;
    parameter Real k3 = 8000.0;
    parameter Real k4 = 40000000.0;
    parameter Real k5 = 1.0;



    Real Br;
    Real BrO3;
    Real Ce;
    Real HBrO2;
    Real HOBr;

initial equation
    Br = 1e-07;
    BrO3 = 0.06;
    Ce = 0.05;
    HBrO2 = 5e-11;
    HOBr = 0.0;

equation

    der(Br) = (f * Ce * k5) - (Br * BrO3 * k1) - (Br * HBrO2 * k2)
    der(BrO3) = 0.0
    der(Ce) = (BrO3 * HBrO2 * k3) - (Ce * k5)
    der(HBrO2) = (Br * BrO3 * k1) + (2.0 * BrO3 * HBrO2 * k3) - (Br * HBrO2 * k2) - (BrO3 * HBrO2 * k3) - (2.0 * pow(HBrO2, 2) * k4)
    der(HOBr) = 0.0

end BIOMD040;
