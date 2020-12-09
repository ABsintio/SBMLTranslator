
model m "Huang1996 - Ultrasensitivity in MAPK cascade"

    parameter Real K_PP_norm_max = 0.900049;
    parameter Real a1 = 1000.0;
    parameter Real d1 = 150.0;
    parameter Real k2 = 150.0;
    parameter Real a2 = 1000.0;
    parameter Real d2 = 150.0;
    parameter Real a3 = 1000.0;
    parameter Real d3 = 150.0;
    parameter Real k3 = 150.0;
    parameter Real a4 = 1000.0;
    parameter Real d4 = 150.0;
    parameter Real k4 = 150.0;
    parameter Real a5 = 1000.0;
    parameter Real d5 = 150.0;
    parameter Real k5 = 150.0;
    parameter Real a6 = 1000.0;
    parameter Real d6 = 150.0;
    parameter Real k6 = 150.0;
    parameter Real a7 = 1000.0;
    parameter Real d7 = 150.0;
    parameter Real k7 = 150.0;
    parameter Real a8 = 1000.0;
    parameter Real d8 = 150.0;
    parameter Real k8 = 150.0;
    parameter Real a9 = 1000.0;
    parameter Real d9 = 150.0;
    parameter Real k9 = 150.0;
    parameter Real a10 = 1000.0;
    parameter Real d10 = 150.0;
    parameter Real k10 = 150.0;



    Real E1;
    Real E2;
    Real KKK;
    Real P_KKK;
    Real KK;
    Real P_KK;
    Real PP_KK;
    Real K;
    Real P_K;
    Real PP_K;
    Real KPase;
    Real KKPase;
    Real E1_KKK;
    Real E2_P_KKK;
    Real P_KKK_KK;
    Real P_KKK_P_KK;
    Real PP_KK_K;
    Real PP_KK_P_K;
    Real KKPase_PP_KK;
    Real KKPase_P_KK;
    Real KPase_PP_K;
    Real KPase_P_K;
    Real K_PP_norm;
    Real KK_PP_norm;
    Real KKK_P_norm;
    Real rel_K_PP_max;

initial equation
    E1 = 3e-05;
    E2 = 0.0003;
    KKK = 0.003;
    P_KKK = 0.0;
    KK = 1.2;
    P_KK = 0.0;
    PP_KK = 0.0;
    K = 1.2;
    P_K = 0.0;
    PP_K = 0.0;
    KPase = 0.12;
    KKPase = 0.0003;
    E1_KKK = 0.0;
    E2_P_KKK = 0.0;
    P_KKK_KK = 0.0;
    P_KKK_P_KK = 0.0;
    PP_KK_K = 0.0;
    PP_KK_P_K = 0.0;
    KKPase_PP_KK = 0.0;
    KKPase_P_KK = 0.0;
    KPase_PP_K = 0.0;
    KPase_P_K = 0.0;
    K_PP_norm = 0.0;
    KK_PP_norm = 0.0;
    KKK_P_norm = 0.0;
    rel_K_PP_max = 0.0;

equation
    K_PP_norm = (PP_K + KPase_PP_K) / (PP_K + P_K + K + PP_KK_K + PP_KK_P_K + KPase_PP_K + KPase_P_K)
    rel_K_PP_max = K_PP_norm / K_PP_norm_max
    KK_PP_norm = (PP_KK + PP_KK_K + PP_KK_P_K + KKPase_PP_KK) / (PP_KK + P_KK + KK + PP_KK_K + PP_KK_P_K + P_KKK_KK + P_KKK_P_KK + KKPase_PP_KK + KKPase_P_KK)
    KKK_P_norm = (P_KKK + P_KKK_KK + P_KKK_P_KK) / (KKK + P_KKK + P_KKK_KK + P_KKK_P_KK)
    der(E1) = (k2 * E1_KKK) - ((a1 * E1 * KKK - d1 * E1_KKK))
    der(E2) = (k2 * E2_P_KKK) - ((a2 * E2 * P_KKK - d2 * E2_P_KKK))
    der(KKK) = (k2 * E2_P_KKK) - ((a1 * E1 * KKK - d1 * E1_KKK))
    der(P_KKK) = (k2 * E1_KKK) + (k3 * P_KKK_KK) + (k5 * P_KKK_P_KK) - ((a2 * E2 * P_KKK - d2 * E2_P_KKK)) - ((a3 * KK * P_KKK - d3 * P_KKK_KK)) - ((a5 * P_KK * P_KKK - d5 * P_KKK_P_KK))
    der(KK) = (k4 * KKPase_P_KK) - ((a3 * KK * P_KKK - d3 * P_KKK_KK))
    der(P_KK) = (k3 * P_KKK_KK) + (k6 * KKPase_PP_KK) - ((a4 * P_KK * KKPase - d4 * KKPase_P_KK)) - ((a5 * P_KK * P_KKK - d5 * P_KKK_P_KK))
    der(PP_KK) = (k5 * P_KKK_P_KK) + (k7 * PP_KK_K) + (k9 * PP_KK_P_K) - ((a6 * PP_KK * KKPase - d6 * KKPase_PP_KK)) - ((a7 * K * PP_KK - d7 * PP_KK_K)) - ((a9 * P_K * PP_KK - d9 * PP_KK_P_K))
    der(K) = (k8 * KPase_P_K) - ((a7 * K * PP_KK - d7 * PP_KK_K))
    der(P_K) = (k7 * PP_KK_K) + (k10 * KPase_PP_K) - ((a8 * P_K * KPase - d8 * KPase_P_K)) - ((a9 * P_K * PP_KK - d9 * PP_KK_P_K))
    der(PP_K) = (k9 * PP_KK_P_K) - ((a10 * PP_K * KPase - d10 * KPase_PP_K))
    der(KPase) = (k8 * KPase_P_K) + (k10 * KPase_PP_K) - ((a8 * P_K * KPase - d8 * KPase_P_K)) - ((a10 * PP_K * KPase - d10 * KPase_PP_K))
    der(KKPase) = (k4 * KKPase_P_KK) + (k6 * KKPase_PP_KK) - ((a4 * P_KK * KKPase - d4 * KKPase_P_KK)) - ((a6 * PP_KK * KKPase - d6 * KKPase_PP_KK))
    der(E1_KKK) = ((a1 * E1 * KKK - d1 * E1_KKK)) - (k2 * E1_KKK)
    der(E2_P_KKK) = ((a2 * E2 * P_KKK - d2 * E2_P_KKK)) - (k2 * E2_P_KKK)
    der(P_KKK_KK) = ((a3 * KK * P_KKK - d3 * P_KKK_KK)) - (k3 * P_KKK_KK)
    der(P_KKK_P_KK) = ((a5 * P_KK * P_KKK - d5 * P_KKK_P_KK)) - (k5 * P_KKK_P_KK)
    der(PP_KK_K) = ((a7 * K * PP_KK - d7 * PP_KK_K)) - (k7 * PP_KK_K)
    der(PP_KK_P_K) = ((a9 * P_K * PP_KK - d9 * PP_KK_P_K)) - (k9 * PP_KK_P_K)
    der(KKPase_PP_KK) = ((a6 * PP_KK * KKPase - d6 * KKPase_PP_KK)) - (k6 * KKPase_PP_KK)
    der(KKPase_P_KK) = ((a4 * P_KK * KKPase - d4 * KKPase_P_KK)) - (k4 * KKPase_P_KK)
    der(KPase_PP_K) = ((a10 * PP_K * KPase - d10 * KPase_PP_K)) - (k10 * KPase_PP_K)
    der(KPase_P_K) = ((a8 * P_K * KPase - d8 * KPase_P_K)) - (k8 * KPase_P_K)

end m;
