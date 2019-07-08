model two_phase_tank

  package Medium = Modelica.Media.R134a.R134a_ph;
  parameter Modelica.SIunits.Volume V = 0.001 "tank volume";
  parameter Modelica.SIunits.Pressure p_inf = 100000 "atm press";
  parameter Modelica.SIunits.Pressure dp_ref = 10000;
  parameter Modelica.SIunits.MassFlowRate mdot_ref = 0.001;
  parameter Modelica.SIunits.HeatFlowRate Q = 10 "heat addition";
  Modelica.SIunits.Mass m "mass";
  Modelica.SIunits.MassFlowRate mdot "mass rate";
  Modelica.SIunits.Pressure p "tank pressure";
  Modelica.SIunits.SpecificEnthalpy h "sp. enthalpy";
  Modelica.SIunits.SpecificEnthalpy u "sp. energy";
  Modelica.SIunits.MassFraction x "quality";
  Modelica.SIunits.Density rho "density";
  Modelica.SIunits.Temperature T "temperature";
  Modelica.SIunits.SpecificEntropy s "entropy";
  
initial equation
  p = 600000;
  h = 300000;

equation
  // dynamics
  mdot = mdot_ref * (p - p_inf) / dp_ref;
  // mass balance
  m = rho * V;
  der(m) = -mdot;
  // energy balance
  der(m*u) = -mdot * h + Q;
  // thermodynamics
  h = u + p / rho;
  x =   Medium.vapourQuality(   Medium.setState_phX(p, h, {1}) );
  rho = Medium.density(         Medium.setState_phX(p, h, {1}) );
  T =   Medium.temperature(     Medium.setState_phX(p, h, {1}) );
  s =   Medium.specificEntropy( Medium.setState_phX(p, h, {1}) );

annotation(
    experiment(StartTime = 0, StopTime = 10, Tolerance = 1e-6, Interval = 0.02));
end two_phase_tank;
