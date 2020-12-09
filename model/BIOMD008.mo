
model m "Gardner1998 - Cell Cycle Goldbeter"

    parameter Real K6 = 0.3;
    parameter Real V1p = 0.75;
    parameter Real V3p = 0.3;
    parameter Real vi = 0.1;
    parameter Real k1 = 0.5;
    parameter Real K5 = 0.02;
    parameter Real kd = 0.02;
    parameter Real K1 = 0.1;
    parameter Real V2 = 0.25;
    parameter Real K2 = 0.1;
    parameter Real K3 = 0.2;
    parameter Real K4 = 0.1;
    parameter Real V4 = 0.1;
    parameter Real a1 = 0.05;
    parameter Real a2 = 0.05;
    parameter Real alpha = 0.1;
    parameter Real d1 = 0.05;
    parameter Real vs = 0.2;

    Real V1(start=0.0);
    Real V3(start=0.0);

    Real C;
    Real X;
    Real M;
    Real Y;
    Real Z;

initial equation
    C = 0.0;
    X = 0.0;
    M = 0.0;
    Y = 0.0;
    Z = 0.0;

equation
    V1 = C * V1p * pow(C + K6, -1)
    V3 = M * V3p
    der(C) = (vi) + (a2 * Z) + (alpha * d1 * Z) - (C * k1 * X * pow(C + K5, -1)) - (C * kd) - (a1 * C * Y)
    der(X) = (V3 * (1 + -1 * X) * pow(K3 + -1 * X + 1, -1)) - (V4 * X * pow(K4 + X, -1))
    der(M) = ((1 + -1 * M) * V1 * pow(K1 + -1 * M + 1, -1)) - (M * V2 * pow(K2 + M, -1))
    der(Y) = (a2 * Z) + (alpha * kd * Z) + (vs) - (a1 * C * Y) - (d1 * Y)
    der(Z) = (a1 * C * Y) - (a2 * Z) - (alpha * d1 * Z) - (alpha * kd * Z)

end m;
