
model BIOMD134 "Izhikevich2004_SpikingNeurons_SpikeLatency"

    function pow
        input  Real x;
        input  Real power;
        output Real y;
        algorithm
            y := x^power;
    end pow;

    parameter Real a = 0.02;
    parameter Real b = 0.2;
    parameter Real c = -65.0;
    parameter Real d = 6.0;
    parameter Real Vthresh = 30.0;

    Real i(start=0.0);
    Real v(start=-70.0);
    Real u(start=-14.0);



initial equation


equation
    when v > Vthresh then
        v = c;
        u = u + d;
    end when;
    when (time > 10) and (time < 13) then
        i = 7.04;
    end when;
    when time >= 13 then
        i = 0;
    end when;


    der(v) = 0.04 * pow(v, 2) + 5 * v + 140 - u + i
    der(u) = a * (b * v - u)

end BIOMD134;
