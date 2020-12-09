
model m "Edelstein1996 - EPSP ACh species"

    parameter Real kf_0 = 300000000.0;
    parameter Real kr_0 = 8000.0;
    parameter Real kf_1 = 150000000.0;
    parameter Real kr_1 = 16000.0;
    parameter Real kf_2 = 30000.0;
    parameter Real kr_2 = 700.0;
    parameter Real kf_3 = 300000000.0;
    parameter Real kr_3 = 8.64;
    parameter Real kf_4 = 150000000.0;
    parameter Real kr_4 = 17.28;
    parameter Real kf_5 = 0.54;
    parameter Real kr_5 = 10800.0;
    parameter Real kf_6 = 130.0;
    parameter Real kr_6 = 2740.0;
    parameter Real kf_7 = 300000000.0;
    parameter Real kr_7 = 4.0;
    parameter Real kf_8 = 150000000.0;
    parameter Real kr_8 = 8.0;
    parameter Real kf_9 = 19.7;
    parameter Real kr_9 = 3.74;
    parameter Real kf_10 = 19.85;
    parameter Real kr_10 = 1.74;
    parameter Real kf_11 = 20.0;
    parameter Real kr_11 = 0.81;
    parameter Real kf_12 = 300000000.0;
    parameter Real kr_12 = 4.0;
    parameter Real kf_13 = 150000000.0;
    parameter Real kr_13 = 8.0;
    parameter Real kf_14 = 0.05;
    parameter Real kr_14 = 0.0012;
    parameter Real kf_15 = 0.05;
    parameter Real kr_15 = 0.0012;
    parameter Real kf_16 = 0.05;
    parameter Real kr_16 = 0.0012;



    Real BLL;
    Real IL;
    Real AL;
    Real A;
    Real BL;
    Real B;
    Real DLL;
    Real D;
    Real ILL;
    Real DL;
    Real I;
    Real ALL;
    Real L;

initial equation
    BLL = 0.0;
    IL = 0.0;
    AL = 0.0;
    A = 0.0;
    BL = 0.0;
    B = 0.0;
    DLL = 0.0;
    D = 0.0;
    ILL = 0.0;
    DL = 0.0;
    I = 0.0;
    ALL = 0.0;
    L = 0.0;

equation

    der(BLL) = ((kf_1 * BL * L - kr_1 * BLL)) - ((kf_2 * BLL - kr_2 * ALL))
    der(IL) = ((kf_7 * I * L - kr_7 * IL)) + ((kf_10 * AL - kr_10 * IL)) - ((kf_8 * IL * L - kr_8 * ILL)) - ((kf_15 * IL - kr_15 * DL))
    der(AL) = ((kf_3 * A * L - kr_3 * AL)) + ((kf_6 * BL - kr_6 * AL)) - ((kf_4 * AL * L - kr_4 * ALL)) - ((kf_10 * AL - kr_10 * IL))
    der(A) = ((kf_5 * B - kr_5 * A)) - ((kf_3 * A * L - kr_3 * AL)) - ((kf_9 * A - kr_9 * I))
    der(BL) = ((kf_0 * B * L - kr_0 * BL)) - ((kf_1 * BL * L - kr_1 * BLL)) - ((kf_6 * BL - kr_6 * AL))
    der(B) =  - ((kf_0 * B * L - kr_0 * BL)) - ((kf_5 * B - kr_5 * A))
    der(DLL) = ((kf_13 * DL * L - kr_13 * DLL)) + ((kf_16 * ILL - kr_16 * DLL)) 
    der(D) = ((kf_14 * I - kr_14 * D)) - ((kf_12 * D * L - kr_12 * DL))
    der(ILL) = ((kf_8 * IL * L - kr_8 * ILL)) + ((kf_11 * ALL - kr_11 * ILL)) - ((kf_16 * ILL - kr_16 * DLL))
    der(DL) = ((kf_12 * D * L - kr_12 * DL)) + ((kf_15 * IL - kr_15 * DL)) - ((kf_13 * DL * L - kr_13 * DLL))
    der(I) = ((kf_9 * A - kr_9 * I)) - ((kf_7 * I * L - kr_7 * IL)) - ((kf_14 * I - kr_14 * D))
    der(ALL) = ((kf_2 * BLL - kr_2 * ALL)) + ((kf_4 * AL * L - kr_4 * ALL)) - ((kf_11 * ALL - kr_11 * ILL))
    der(L) =  - ((kf_0 * B * L - kr_0 * BL)) - ((kf_1 * BL * L - kr_1 * BLL)) - ((kf_3 * A * L - kr_3 * AL)) - ((kf_4 * AL * L - kr_4 * ALL)) - ((kf_7 * I * L - kr_7 * IL)) - ((kf_8 * IL * L - kr_8 * ILL)) - ((kf_12 * D * L - kr_12 * DL)) - ((kf_13 * DL * L - kr_13 * DLL))

end m;
