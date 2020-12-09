model BIOMD069 "Fuss2006_MitoticActivation"

    parameter Real k1 = 1.0;
    parameter Real k2 = 0.8;
    parameter Real k3 = 1.0;
    parameter Real k4 = 10.0;
    parameter Real kPTP = 1.0;
    parameter Real kCbp = 1.0;
    parameter Real p1 = 0.05;
    parameter Real p2 = 0.15;
    parameter Real p3 = 0.035;
    parameter Real src_background = 0.0001;
    parameter Real PTP_background = 0.0;
    parameter Real kCSKon = 0.1;
    parameter Real kCSKoff = 0.01;
    parameter Real rho_srca = 1.0;
    parameter Real rho_srco = 0.0;
    parameter Real rho_srcc = 1.0;
    parameter Real Kser = 1.0;
    parameter Real acsk0 = 0.0;

    Real ptp_activity(start=0.0);
    Real src_activity(start=0.0);

    Real srci;
    Real srco;
    Real srca;
    Real srcc;
    Real Cbp_P_CSK;
    Real CSK_cytoplasm;
    Real PTP;
    Real PTP_pY789;
    Real Cbp;
    Real Cbp_P;

initial equation
    srci = 1.0;
    srco = 0.0;
    srca = 0.0;
    srcc = 0.0;
    Cbp_P_CSK = 0.0;
    CSK_cytoplasm = 1.0;
    PTP = 1.0;
    PTP_pY789 = 0.0;
    Cbp = 1.0;
    Cbp_P = 0.0;

equation
    src_activity = rho_srco * srco + rho_srca * srca + src_background + rho_srcc * srcc
    ptp_activity = PTP_background + Kser * PTP_pY789
    der(srci) = (k4 * p1 * srcc) - ((k2 * ptp_activity * srci - k1 * Cbp_P_CSK * srco))
    der(srco) = ((k2 * ptp_activity * srci - k1 * Cbp_P_CSK * srco)) - ((k3 * src_activity * srco - p1 * srca))
    der(srca) = ((k3 * src_activity * srco - p1 * srca)) - ((k1 * Cbp_P_CSK * srca - k2 * ptp_activity * srcc))
    der(srcc) = ((k1 * Cbp_P_CSK * srca - k2 * ptp_activity * srcc)) - (k4 * p1 * srcc)
    der(Cbp_P_CSK) = ((Cbp_P * kCSKon * CSK_cytoplasm - kCSKoff * Cbp_P_CSK)) 
    der(CSK_cytoplasm) =  - ((Cbp_P * kCSKon * CSK_cytoplasm - kCSKoff * Cbp_P_CSK))
    der(PTP) =  - (((kPTP * src_activity + p3) * PTP - p2 * PTP_pY789))
    der(PTP_pY789) = (((kPTP * src_activity + p3) * PTP - p2 * PTP_pY789)) 
    der(Cbp) =  - (kCbp * src_activity * Cbp)
    der(Cbp_P) = (kCbp * src_activity * Cbp) - ((Cbp_P * kCSKon * CSK_cytoplasm - kCSKoff * Cbp_P_CSK))

end BIOMD069;
