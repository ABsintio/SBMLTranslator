
model m "Tyson1991 - Cell Cycle 6 var"

    parameter Real k6 = 1.0;
    parameter Real k8notP = 1000000.0;
    parameter Real k9 = 1000.0;
    parameter Real k3 = 200.0;
    parameter Real k5notP = 0.0;
    parameter Real k1aa = 0.015;
    parameter Real k2 = 0.0;
    parameter Real k7 = 0.6;
    parameter Real k4 = 180.0;
    parameter Real k4prime = 0.018;



    Real EmptySet;
    Real C2;
    Real CP;
    Real M;
    Real pM;
    Real Y;
    Real YP;
    Real YT;
    Real CT;

initial equation
    EmptySet = 0.0;
    C2 = 0.0;
    CP = 0.0;
    M = 0.0;
    pM = 0.0;
    Y = 0.0;
    YP = 0.0;
    YT = 0.0;
    CT = 0.0;

equation
    YT = Y + YP + M + pM
    CT = C2 + CP + M + pM
    der(EmptySet) = 0.0
    der(C2) = (k6 * M) + (CP * k9) - (C2 * k8notP)
    der(CP) = (C2 * k8notP) - (CP * k9) - (CP * k3 * Y)
    der(M) = (pM * (k4prime + k4 * pow(M / CT, 2))) - (k6 * M) - (k5notP * M)
    der(pM) = (CP * k3 * Y) + (k5notP * M) - (pM * (k4prime + k4 * pow(M / CT, 2)))
    der(Y) = (k1aa) - (CP * k3 * Y) - (k2 * Y)
    der(YP) = (k6 * M) - (k7 * YP)

end m;
