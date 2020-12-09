
model m "Kholodenko2000 - Ultrasensitivity and negative feedback bring oscillations in MAPK cascade"

    parameter Real V1 = 2.5;
    parameter Real Ki = 9.0;
    parameter Real n = 1.0;
    parameter Real K1 = 10.0;
    parameter Real V2 = 0.25;
    parameter Real KK2 = 8.0;
    parameter Real k3 = 0.025;
    parameter Real KK3 = 15.0;
    parameter Real k4 = 0.025;
    parameter Real KK4 = 15.0;
    parameter Real V5 = 0.75;
    parameter Real KK5 = 15.0;
    parameter Real V6 = 0.75;
    parameter Real KK6 = 15.0;
    parameter Real k7 = 0.025;
    parameter Real KK7 = 15.0;
    parameter Real k8 = 0.025;
    parameter Real KK8 = 15.0;
    parameter Real V9 = 0.5;
    parameter Real KK9 = 15.0;
    parameter Real V10 = 0.5;
    parameter Real KK10 = 15.0;



    Real MKKK;
    Real MKKK_P;
    Real MKK;
    Real MKK_P;
    Real MKK_PP;
    Real MAPK;
    Real MAPK_P;
    Real MAPK_PP;

initial equation
    MKKK = 90.0;
    MKKK_P = 10.0;
    MKK = 280.0;
    MKK_P = 10.0;
    MKK_PP = 10.0;
    MAPK = 280.0;
    MAPK_P = 10.0;
    MAPK_PP = 10.0;

equation

    der(MKKK) = (V2 * MKKK_P / (KK2 + MKKK_P)) - (V1 * MKKK / ((1 + pow(MAPK_PP / Ki, n)) * (K1 + MKKK)))
    der(MKKK_P) = (V1 * MKKK / ((1 + pow(MAPK_PP / Ki, n)) * (K1 + MKKK))) - (V2 * MKKK_P / (KK2 + MKKK_P))
    der(MKK) = (V6 * MKK_P / (KK6 + MKK_P)) - (k3 * MKKK_P * MKK / (KK3 + MKK))
    der(MKK_P) = (k3 * MKKK_P * MKK / (KK3 + MKK)) + (V5 * MKK_PP / (KK5 + MKK_PP)) - (k4 * MKKK_P * MKK_P / (KK4 + MKK_P)) - (V6 * MKK_P / (KK6 + MKK_P))
    der(MKK_PP) = (k4 * MKKK_P * MKK_P / (KK4 + MKK_P)) - (V5 * MKK_PP / (KK5 + MKK_PP))
    der(MAPK) = (V10 * MAPK_P / (KK10 + MAPK_P)) - (k7 * MKK_PP * MAPK / (KK7 + MAPK))
    der(MAPK_P) = (k7 * MKK_PP * MAPK / (KK7 + MAPK)) + (V9 * MAPK_PP / (KK9 + MAPK_PP)) - (k8 * MKK_PP * MAPK_P / (KK8 + MAPK_P)) - (V10 * MAPK_P / (KK10 + MAPK_P))
    der(MAPK_PP) = (k8 * MKK_PP * MAPK_P / (KK8 + MAPK_P)) - (V9 * MAPK_PP / (KK9 + MAPK_PP))

end m;
