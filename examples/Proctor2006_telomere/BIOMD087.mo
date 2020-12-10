
model BIOMD087 "Proctor2006_telomere"

    parameter Real k1 = 0.0005;
    parameter Real k2 = 0.000385;
    parameter Real k3 = 1.5e-08;
    parameter Real k4 = 0.01;
    parameter Real k5 = 0.0003;
    parameter Real k6a = 5e-05;
    parameter Real k6b = 0.0005;
    parameter Real k7a = 3e-05;
    parameter Real k7b = 3e-05;
    parameter Real k8a = 0.001;
    parameter Real k8b = 100.0;
    parameter Real k8c = 100.0;
    parameter Real k8d = 0.004;
    parameter Real k9 = 100.0;
    parameter Real k10a = 0.05;
    parameter Real k10b = 0.05;
    parameter Real k11 = 1e-05;
    parameter Real k12 = 0.00017;
    parameter Real k13 = 1.0;
    parameter Real k14 = 3.3e-06;
    parameter Real k15 = 0.2;
    parameter Real k16 = 0.1;
    parameter Real k17a = 0.05;
    parameter Real k17b = 0.05;
    parameter Real k18a = 0.001;
    parameter Real k18b = 1e-05;
    parameter Real k19 = 0.001;
    parameter Real kc1 = 0.16;
    parameter Real kc2 = 0.01;
    parameter Real kc3 = 0.0012;
    parameter Real kc4 = 0.01;

    Real kalive(start=1.0);

    Real Ctelo;
    Real Utelo;
    Real Cdc13;
    Real Rad17Utelo;
    Real Rad17;
    Real Rad24;
    Real RPA;
    Real Mec1;
    Real ssDNA;
    Real RPAssDNA;
    Real RPAssDNA1;
    Real RPAssDNA2;
    Real Mec1RPAssDNA;
    Real ExoXI;
    Real ExoXA;
    Real Exo1I;
    Real Exo1A;
    Real Rad9I;
    Real Rad9A;
    Real Rad53I;
    Real Rad53A;
    Real Chk1I;
    Real Chk1A;
    Real Dun1I;
    Real Dun1A;
    Real ATP;
    Real ADP;
    Real Rad9Kin;
    Real recovery;
    Real G1;
    Real S;
    Real G2;
    Real M;
    Real G1cyclin;
    Real Scyclin;
    Real G2cyclin;
    Real Mcyclin;
    Real G1CdkI;
    Real G1CdkA;
    Real SCdkI;
    Real SCdkA;
    Real G2CdkI;
    Real G2CdkA;
    Real MCdkI;
    Real MCdkA;
    Real G1Soff;
    Real G1Son;
    Real SG2off;
    Real SG2on;
    Real G2Moff;
    Real G2Mon;
    Real MG1off;
    Real MG1on;
    Real sink;
    Real budscar;

initial equation
    Ctelo = 0.0;
    Utelo = 0.0;
    Cdc13 = 0.0;
    Rad17Utelo = 0.0;
    Rad17 = 0.0;
    Rad24 = 0.0;
    RPA = 0.0;
    Mec1 = 0.0;
    ssDNA = 0.0;
    RPAssDNA = 0.0;
    RPAssDNA1 = 0.0;
    RPAssDNA2 = 0.0;
    Mec1RPAssDNA = 0.0;
    ExoXI = 0.0;
    ExoXA = 0.0;
    Exo1I = 0.0;
    Exo1A = 0.0;
    Rad9I = 0.0;
    Rad9A = 0.0;
    Rad53I = 0.0;
    Rad53A = 0.0;
    Chk1I = 0.0;
    Chk1A = 0.0;
    Dun1I = 0.0;
    Dun1A = 0.0;
    ATP = 0.0;
    ADP = 0.0;
    Rad9Kin = 0.0;
    recovery = 0.0;
    G1 = 0.0;
    S = 0.0;
    G2 = 0.0;
    M = 0.0;
    G1cyclin = 0.0;
    Scyclin = 0.0;
    G2cyclin = 0.0;
    Mcyclin = 0.0;
    G1CdkI = 0.0;
    G1CdkA = 0.0;
    SCdkI = 0.0;
    SCdkA = 0.0;
    G2CdkI = 0.0;
    G2CdkA = 0.0;
    MCdkI = 0.0;
    MCdkA = 0.0;
    G1Soff = 0.0;
    G1Son = 0.0;
    SG2off = 0.0;
    SG2on = 0.0;
    G2Moff = 0.0;
    G2Mon = 0.0;
    MG1off = 0.0;
    MG1on = 0.0;
    sink = 0.0;
    budscar = 0.0;

equation
    when Mec1RPAssDNA >= 800 then
        Rad9Kin = 1;
    end when;
    when (Mec1RPAssDNA + RPAssDNA + ssDNA) <= 1 then
        recovery = 1;
        Mec1RPAssDNA = 0;
        RPAssDNA = 0;
        ssDNA = 0;
    end when;
    when (G2 == 1) and (Rad17Utelo == 0) then
        G2Mon = 1;
        G2Moff = 0;
        recovery = 0;
        Rad9A = 0;
        Rad9I = 20;
        Chk1A = 0;
        Chk1I = 60;
        Dun1A = 0;
        Dun1I = 3000;
        Exo1A = 0;
        Exo1I = 670;
        ExoXA = 0;
        ExoXI = 70;
        Rad53I = 6900;
        Rad53A = 0;
    end when;
    when Rad17Utelo == 0 then
        recovery = 0;
    end when;
    when G1cyclin > 100 then
        G1CdkA = 1;
        G1CdkI = 0;
    end when;
    when Scyclin > 100 then
        SCdkA = 1;
        SCdkI = 0;
    end when;
    when G2cyclin > 100 then
        G2CdkA = 1;
        G2CdkI = 0;
    end when;
    when Mcyclin > 100 then
        MCdkA = 1;
        MCdkI = 0;
    end when;
    when (Mec1RPAssDNA + RPAssDNA + ssDNA) >= 2000 then
        kalive = 0;
    end when;


    der(Ctelo) = (k1 * Cdc13 * Utelo * kalive) + (Cdc13 * k19 * Rad17Utelo * recovery * kalive) - (k2 * Ctelo * kalive)
    der(Utelo) = (k2 * Ctelo * kalive) + (k7a * Utelo * Exo1A * kalive) - (k1 * Cdc13 * Utelo * kalive) - (k3 * Utelo * Rad17 * Rad24 * ATP * kalive / (5000 + ATP)) - (k7a * Utelo * Exo1A * kalive)
    der(Cdc13) = (k2 * Ctelo * kalive) - (k1 * Cdc13 * Utelo * kalive) - (Cdc13 * k19 * Rad17Utelo * recovery * kalive)
    der(Rad17Utelo) = (k3 * Utelo * Rad17 * Rad24 * ATP * kalive / (5000 + ATP)) + (k4 * ExoXI * Rad17Utelo * kalive) + (k5 * ExoXA * Rad17Utelo * kalive) + (k7b * Rad17Utelo * Exo1A * kalive) - (k4 * ExoXI * Rad17Utelo * kalive) - (k5 * ExoXA * Rad17Utelo * kalive) - (k7b * Rad17Utelo * Exo1A * kalive) - (Cdc13 * k19 * Rad17Utelo * recovery * kalive)
    der(Rad17) = (Cdc13 * k19 * Rad17Utelo * recovery * kalive) - (k3 * Utelo * Rad17 * Rad24 * ATP * kalive / (5000 + ATP))
    der(Rad24) = (k3 * Utelo * Rad17 * Rad24 * ATP * kalive / (5000 + ATP)) + (k6b * Exo1I * Rad24 * kalive) - (k3 * Utelo * Rad17 * Rad24 * ATP * kalive / (5000 + ATP)) - (k6b * Exo1I * Rad24 * kalive)
    der(RPA) = (k17a * Mec1RPAssDNA * S * kalive) + (G2 * G2Moff * k17b * Mec1RPAssDNA * kalive) - (k8a * RPA * ssDNA * kalive)
    der(Mec1) = (k17a * Mec1RPAssDNA * S * kalive) + (G2 * G2Moff * k17b * Mec1RPAssDNA * kalive) - (k8d * RPAssDNA * Mec1 * kalive)
    der(ssDNA) = (k5 * ExoXA * Rad17Utelo * kalive) + (k7a * Utelo * Exo1A * kalive) + (k7b * Rad17Utelo * Exo1A * kalive) + (3.0 * k17a * Mec1RPAssDNA * S * kalive) + (3.0 * G2 * G2Moff * k17b * Mec1RPAssDNA * kalive) - (k8a * RPA * ssDNA * kalive) - (k8b * RPAssDNA1 * ssDNA * kalive) - (k8c * RPAssDNA2 * ssDNA * kalive) - (k18a * S * ssDNA * kalive) - (G2 * G2Moff * k18b * ssDNA * kalive)
    der(RPAssDNA) = (k8c * RPAssDNA2 * ssDNA * kalive) - (k8d * RPAssDNA * Mec1 * kalive)
    der(RPAssDNA1) = (k8a * RPA * ssDNA * kalive) - (k8b * RPAssDNA1 * ssDNA * kalive)
    der(RPAssDNA2) = (k8b * RPAssDNA1 * ssDNA * kalive) - (k8c * RPAssDNA2 * ssDNA * kalive)
    der(Mec1RPAssDNA) = (k8d * RPAssDNA * Mec1 * kalive) - (k17a * Mec1RPAssDNA * S * kalive) - (G2 * G2Moff * k17b * Mec1RPAssDNA * kalive)
    der(ExoXI) = (ExoXA * k10a * Rad9A * kalive) + (ExoXA * k10b * Rad9I * kalive) - (k4 * ExoXI * Rad17Utelo * kalive)
    der(ExoXA) = (k4 * ExoXI * Rad17Utelo * kalive) + (k5 * ExoXA * Rad17Utelo * kalive) - (k5 * ExoXA * Rad17Utelo * kalive) - (ExoXA * k10a * Rad9A * kalive) - (ExoXA * k10b * Rad9I * kalive)
    der(Exo1I) = (Exo1A * k13 * Rad53A * kalive) - (k6a * Exo1I * kalive) - (k6b * Exo1I * Rad24 * kalive)
    der(Exo1A) = (k6a * Exo1I * kalive) + (k6b * Exo1I * Rad24 * kalive) + (k7a * Utelo * Exo1A * kalive) + (k7b * Rad17Utelo * Exo1A * kalive) - (k7a * Utelo * Exo1A * kalive) - (k7b * Rad17Utelo * Exo1A * kalive) - (Exo1A * k13 * Rad53A * kalive)
    der(Rad9I) = (ExoXA * k10b * Rad9I * kalive) - (k9 * Rad9Kin * Rad9I * kalive) - (ExoXA * k10b * Rad9I * kalive)
    der(Rad9A) = (k9 * Rad9Kin * Rad9I * kalive) + (ExoXA * k10a * Rad9A * kalive) + (k11 * Rad53I * Rad9A * kalive) + (Chk1I * k12 * Rad9A * kalive) - (ExoXA * k10a * Rad9A * kalive) - (k11 * Rad53I * Rad9A * kalive) - (Chk1I * k12 * Rad9A * kalive)
    der(Rad53I) =  - (k11 * Rad53I * Rad9A * kalive)
    der(Rad53A) = (k11 * Rad53I * Rad9A * kalive) + (Exo1A * k13 * Rad53A * kalive) + (Dun1I * k14 * Rad53A * kalive) - (Exo1A * k13 * Rad53A * kalive) - (Dun1I * k14 * Rad53A * kalive)
    der(Chk1I) =  - (Chk1I * k12 * Rad9A * kalive)
    der(Chk1A) = (Chk1I * k12 * Rad9A * kalive) + (Chk1A * G2Mon * k15 * kalive) - (Chk1A * G2Mon * k15 * kalive)
    der(Dun1I) =  - (Dun1I * k14 * Rad53A * kalive)
    der(Dun1A) = (Dun1I * k14 * Rad53A * kalive) + (Dun1A * G2Mon * k16 * kalive) - (Dun1A * G2Mon * k16 * kalive)
    der(ATP) =  - (k3 * Utelo * Rad17 * Rad24 * ATP * kalive / (5000 + ATP))
    der(ADP) = (k3 * Utelo * Rad17 * Rad24 * ATP * kalive / (5000 + ATP)) 
    der(Rad9Kin) = (k9 * Rad9Kin * Rad9I * kalive) - (k9 * Rad9Kin * Rad9I * kalive)
    der(recovery) = (Cdc13 * k19 * Rad17Utelo * recovery * kalive) - (Cdc13 * k19 * Rad17Utelo * recovery * kalive)
    der(G1) = (G1 * kc1 * kalive) + (G1 * G1CdkA * G1Soff * kc2 * kalive) + (kc4 * M * MCdkA * MG1on * kalive) - (G1 * kc1 * kalive) - (G1 * G1CdkA * G1Soff * kc2 * kalive) - (G1 * G1CdkA * G1Son * kc4 * kalive)
    der(S) = (k17a * Mec1RPAssDNA * S * kalive) + (k18a * S * ssDNA * kalive) + (kc1 * S * kalive) + (kc2 * S * SCdkA * SG2off * kalive) + (G1 * G1CdkA * G1Son * kc4 * kalive) - (k17a * Mec1RPAssDNA * S * kalive) - (k18a * S * ssDNA * kalive) - (kc1 * S * kalive) - (kc2 * S * SCdkA * SG2off * kalive) - (kc4 * S * SCdkA * SG2on * kalive)
    der(G2) = (G2 * G2Moff * k17b * Mec1RPAssDNA * kalive) + (G2 * G2Moff * k18b * ssDNA * kalive) + (G2 * kc1 * kalive) + (G2 * G2CdkA * G2Moff * kc2 * kalive) + (kc4 * S * SCdkA * SG2on * kalive) - (G2 * G2Moff * k17b * Mec1RPAssDNA * kalive) - (G2 * G2Moff * k18b * ssDNA * kalive) - (G2 * kc1 * kalive) - (G2 * G2CdkA * G2Moff * kc2 * kalive) - (G2 * G2CdkA * G2Mon * kc4 * kalive)
    der(M) = (kc1 * M * kalive) + (kc2 * M * MCdkA * MG1off * kalive) + (G2 * G2CdkA * G2Mon * kc4 * kalive) - (kc1 * M * kalive) - (kc2 * M * MCdkA * MG1off * kalive) - (kc4 * M * MCdkA * MG1on * kalive)
    der(G1cyclin) = (G1 * kc1 * kalive) - (G1cyclin * kc3 * kalive)
    der(Scyclin) = (kc1 * S * kalive) - (kc3 * Scyclin * kalive)
    der(G2cyclin) = (G2 * kc1 * kalive) - (G2cyclin * kc3 * kalive)
    der(Mcyclin) = (kc1 * M * kalive) - (kc3 * Mcyclin * kalive)
    der(G1CdkI) = (G1 * G1CdkA * G1Son * kc4 * kalive) 
    der(G1CdkA) = (G1 * G1CdkA * G1Soff * kc2 * kalive) - (G1 * G1CdkA * G1Soff * kc2 * kalive) - (G1 * G1CdkA * G1Son * kc4 * kalive)
    der(SCdkI) = (kc4 * S * SCdkA * SG2on * kalive) 
    der(SCdkA) = (kc2 * S * SCdkA * SG2off * kalive) - (kc2 * S * SCdkA * SG2off * kalive) - (kc4 * S * SCdkA * SG2on * kalive)
    der(G2CdkI) = (G2 * G2CdkA * G2Mon * kc4 * kalive) 
    der(G2CdkA) = (G2 * G2CdkA * G2Moff * kc2 * kalive) - (G2 * G2CdkA * G2Moff * kc2 * kalive) - (G2 * G2CdkA * G2Mon * kc4 * kalive)
    der(MCdkI) = (kc4 * M * MCdkA * MG1on * kalive) 
    der(MCdkA) = (kc2 * M * MCdkA * MG1off * kalive) - (kc2 * M * MCdkA * MG1off * kalive) - (kc4 * M * MCdkA * MG1on * kalive)
    der(G1Soff) = (G1 * G1CdkA * G1Son * kc4 * kalive) - (G1 * G1CdkA * G1Soff * kc2 * kalive)
    der(G1Son) = (G1 * G1CdkA * G1Soff * kc2 * kalive) - (G1 * G1CdkA * G1Son * kc4 * kalive)
    der(SG2off) = (kc4 * S * SCdkA * SG2on * kalive) - (kc2 * S * SCdkA * SG2off * kalive)
    der(SG2on) = (kc2 * S * SCdkA * SG2off * kalive) - (kc4 * S * SCdkA * SG2on * kalive)
    der(G2Moff) = (Chk1A * G2Mon * k15 * kalive) + (Dun1A * G2Mon * k16 * kalive) + (G2 * G2Moff * k17b * Mec1RPAssDNA * kalive) + (G2 * G2Moff * k18b * ssDNA * kalive) + (G2 * G2CdkA * G2Mon * kc4 * kalive) - (G2 * G2Moff * k17b * Mec1RPAssDNA * kalive) - (G2 * G2Moff * k18b * ssDNA * kalive) - (G2 * G2CdkA * G2Moff * kc2 * kalive)
    der(G2Mon) = (G2 * G2CdkA * G2Moff * kc2 * kalive) - (Chk1A * G2Mon * k15 * kalive) - (Dun1A * G2Mon * k16 * kalive) - (G2 * G2CdkA * G2Mon * kc4 * kalive)
    der(MG1off) = (kc4 * M * MCdkA * MG1on * kalive) - (kc2 * M * MCdkA * MG1off * kalive)
    der(MG1on) = (kc2 * M * MCdkA * MG1off * kalive) - (kc4 * M * MCdkA * MG1on * kalive)
    der(sink) = 0.0
    der(budscar) = (kc4 * M * MCdkA * MG1on * kalive) 

end BIOMD087;
