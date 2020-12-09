
model m "Novak1997 - Cell Cycle"

    parameter Real mu = 0.00495;
    parameter Real k1 = 0.015;
    parameter Real k2prime = 0.05;
    parameter Real k3 = 0.09375;
    parameter Real k4 = 0.1875;
    parameter Real k5 = 0.00175;
    parameter Real k6prime = 0.0;
    parameter Real k7 = 100.0;
    parameter Real k7r = 0.1;
    parameter Real k8 = 10.0;
    parameter Real k8r = 0.1;
    parameter Real kc = 1.0;
    parameter Real kcr = 0.25;
    parameter Real ki = 0.4;
    parameter Real kir = 0.1;
    parameter Real ku = 0.2;
    parameter Real kur = 0.1;
    parameter Real ku2 = 1.0;
    parameter Real kur2 = 0.3;
    parameter Real kw = 1.0;
    parameter Real kwr = 0.25;
    parameter Real V2 = 0.25;
    parameter Real V2prime = 0.0075;
    parameter Real V6 = 7.5;
    parameter Real V6prime = 0.0375;
    parameter Real V25 = 0.5;
    parameter Real V25prime = 0.025;
    parameter Real Vw = 0.35;
    parameter Real Vwprime = 0.035;
    parameter Real Kmc = 0.1;
    parameter Real Kmcr = 0.1;
    parameter Real Kmi = 0.01;
    parameter Real Kmir = 0.01;
    parameter Real Kmp = 0.001;
    parameter Real Kmu = 0.01;
    parameter Real Kmur = 0.01;
    parameter Real Kmu2 = 0.05;
    parameter Real Kmur2 = 0.05;
    parameter Real Kmw = 0.1;
    parameter Real Kmwr = 0.1;
    parameter Real alpha = 0.25;
    parameter Real beta = 0.05;
    parameter Real Cig1 = 0.0;

    Real Mass(start=0.49);
    Real kp(start=3.25);
    Real k2(start=0.0);
    Real k6(start=0.0);
    Real kwee(start=0.0);
    Real k25(start=0.0);

    Real UbE;
    Real UbE2;
    Real Wee1;
    Real Cdc25;
    Real G2K;
    Real R;
    Real G1K;
    Real IE;
    Real PG2;
    Real G1R;
    Real G2R;
    Real PG2R;
    Real SPF;
    Real MPF;
    Real IEB;
    Real UbEB;
    Real UbE2B;
    Real Wee1B;
    Real Cdc25B;
    Real Rum1Total;
    Real Cdc13Total;
    Real Cig2Total;

initial equation
    UbE = 0.0;
    UbE2 = 0.0;
    Wee1 = 0.0;
    Cdc25 = 0.0;
    G2K = 0.0;
    R = 0.0;
    G1K = 0.0;
    IE = 0.0;
    PG2 = 0.0;
    G1R = 0.0;
    G2R = 0.0;
    PG2R = 0.0;
    SPF = 0.0;
    MPF = 0.0;
    IEB = 0.0;
    UbEB = 0.0;
    UbE2B = 0.0;
    Wee1B = 0.0;
    Cdc25B = 0.0;
    Rum1Total = 0.0;
    Cdc13Total = 0.0;
    Cig2Total = 0.0;

equation
    IEB = 1 - IE
    UbEB = 1 - UbE
    UbE2B = 1 - UbE2
    Wee1B = 1 - Wee1
    Cdc25B = 1 - Cdc25
    Rum1Total = G1R + G2R + PG2R + R
    Cdc13Total = G2K + G2R + PG2 + PG2R
    Cig2Total = G1K + G1R
    k2 = UbE * V2 + (1 - UbE) * V2prime
    k6 = UbE2 * V6 + (1 - UbE2) * V6prime
    kwee = Vwprime * (1 - Wee1) + Vw * Wee1
    k25 = Cdc25 * V25 + (1 - Cdc25) * V25prime
    MPF = G2K + beta * PG2
    SPF = Cig1 + alpha * G1K + MPF
    der(Mass) = Mass * mu
    der(UbE) = (IE * ku * UbEB / (Kmu + UbEB) - kur * UbE / (Kmur + UbE)) 
    der(UbE2) = (ku2 * MPF * UbE2B / (Kmu2 + UbE2B) - kur2 * UbE2 / (Kmur2 + UbE2)) 
    der(Wee1) = (kwr * Wee1B / (Kmwr + Wee1B) - kw * MPF * Wee1 / (Kmw + Wee1)) 
    der(Cdc25) = (Cdc25B * kc * MPF / (Cdc25B + Kmc) - Cdc25 * kcr / (Cdc25 + Kmcr)) 
    der(G2K) = (k1) + (G2R * k4) - (G2K * kwee - k25 * PG2) - (G2K * k7 * R - G2R * k7r) - (G2K * k2)
    der(R) = (k2prime * PG2R) + (G2R * k2prime) + (G1R * k6prime) + (G2R * k2) + (k2 * PG2R) + (k3) - (G2K * k7 * R - G2R * k7r) - (k7 * PG2 * R - k7r * PG2R) - (k4 * R) - (G1K * k8 * R - G1R * k8r) - (kp * Mass * R * SPF / (Kmp + R))
    der(G1K) = (k5) + (G1R * k4) - (G1K * k6) - (G1K * k8 * R - G1R * k8r)
    der(IE) = (IEB * ki * MPF / (IEB + Kmi) - IE * kir / (IE + Kmir)) 
    der(PG2) = (G2K * kwee - k25 * PG2) + (k4 * PG2R) - (k7 * PG2 * R - k7r * PG2R) - (k2 * PG2)
    der(G1R) = (G1K * k8 * R - G1R * k8r) - (G1R * k4) - (G1R * k6prime)
    der(G2R) = (G2K * k7 * R - G2R * k7r) - (G2R * k4) - (G2R * k2prime) - (G2R * k2)
    der(PG2R) = (k7 * PG2 * R - k7r * PG2R) - (k4 * PG2R) - (k2prime * PG2R) - (k2 * PG2R)

end m;
