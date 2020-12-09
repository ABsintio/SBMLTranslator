
model m "Goldbeter1991 - Min Mit Oscil, Expl Inact"

    parameter Real VM1 = 3.0;
    parameter Real VM3 = 1.0;
    parameter Real Kc = 0.5;
    parameter Real vi = 0.025;
    parameter Real kd = 0.01;
    parameter Real vd = 0.25;
    parameter Real Kd = 0.02;
    parameter Real K1 = 0.005;
    parameter Real V2 = 1.5;
    parameter Real K2 = 0.005;
    parameter Real K3 = 0.005;
    parameter Real K4 = 0.005;
    parameter Real V4 = 0.5;

    Real V1(start=0.0);
    Real V3(start=0.0);

    Real C;
    Real M;
    Real X;
    Real MI;
    Real XI;

initial equation
    C = 0.01;
    M = 0.01;
    X = 0.01;
    MI = 0.99;
    XI = 0.99;

equation
    V1 = C * VM1 * pow(C + Kc, -1)
    V3 = M * VM3
    der(C) = (vi) - (C * kd) - (C * vd * X * pow(C + Kd, -1))
    der(M) = (MI * V1 * pow(K1 + MI, -1)) - (M * V2 * pow(K2 + M, -1))
    der(X) = (V3 * XI * pow(K3 + XI, -1)) - (V4 * X * pow(K4 + X, -1))
    der(MI) = (M * V2 * pow(K2 + M, -1)) - (MI * V1 * pow(K1 + MI, -1))
    der(XI) = (V4 * X * pow(K4 + X, -1)) - (V3 * XI * pow(K3 + XI, -1))

end m;
