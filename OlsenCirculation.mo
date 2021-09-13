within ;
package OlsenCirculation "A reimplementation in Modelica language"
  partial model partialDriving_Olsen
    "Model for atrial and ventricular activation"
    parameter Boolean fixedActivation = false;
    parameter Real fixedActValue = 0;
  // Heart rate
    parameter Boolean HR_stepUp=false;
    parameter Real stepUpLength = 10;
    parameter Physiolibrary.Types.Frequency HRo0_ref(displayUnit="1/min")=
      1.0333333333333
      "Reference cycle length for scaling the timing";
    parameter Physiolibrary.Types.Frequency HRo0(displayUnit="1/min")=
      1.0333333333333;
    Physiolibrary.Types.Frequency HRo "Olsen's heart rate";
    Physiolibrary.Types.Fraction cycleScale=HRo0_ref/HRo;
  // Activation, global
  parameter Modelica.Units.SI.Time Atact(displayUnit="ms")=0.01
                                              " Timing of atrial activation within cycle";
  parameter Modelica.Units.SI.Time Aup(displayUnit="ms")=0.1
                                             " Duration of upstroke of atrial activation";
  parameter Modelica.Units.SI.Time Aactdur(displayUnit="ms")=0.3
                                                 " Duration of atrial activation";
  parameter Modelica.Units.SI.Time Adown(displayUnit="ms")=0.16
                                               " Duration of downstroke of atrial activation";

  parameter Real PRRRrel = 0.05 " Regression of PR change vs. RR change";
   // Ref: Subject-specific heart rate dependency of electrocardiographic QT, PQ, and QRS intervals.
   // Malik M1, Hnatkova K, Sisakova M, Schmidt G.
  parameter Modelica.Units.SI.Time LVdelay0(displayUnit="ms")=0.151
                                                  " Delay btw. atria and LV at HR 60";
    Real LVdelay=LVdelay0 + PRRRrel*((1/HRo) - 1)
      " Delay btw atria and LV with HR dependence";
  Real LVtact = Atact + LVdelay " Timing of LV activation within cycle";
  parameter Modelica.Units.SI.Time LVup(displayUnit="ms")=0.243
                                              " Duration of upstroke of LV activation - no relation with HR (J Electrocardiol. 2008 Nov-Dec;41(6):491-7. doi: 10.1016/j.jelectrocard.2008.06.022. Epub 2008 Sep 24.";
  parameter Modelica.Units.SI.Time LVactdur0(displayUnit="ms")=0.571
                                                   " Duration of LV activation at HR 60";
  parameter Real LVdownpart=0.5                      " Duration of downstroke of LV activation as fraction of total duration";
    Real LVactdur=LVactdur0*(((1)/HRo)^(1/3))
      " Duration of LV activation, HR dependent";
  // Ref: Fridericia LS (1920). "The duration of systole in the electrocardiogram of normal subjects and of patients with heart disease". Acta Medica Scandinavica (53): 469â€“486.
  Real LVdown = LVactdur*LVdownpart " Duration of downstroke of LV activation";

  parameter Modelica.Units.SI.Time RVdelay0(displayUnit="ms")=0.181
                                                  " Delay btw. atria and RV at HR 60";
    Real RVdelay=RVdelay0 + PRRRrel*((1/HRo) - 1)
      " Delay btw atria and LV, HR dependent";
  Real RVtact = Atact + RVdelay " Timing of RV activation within cycle";
  parameter Modelica.Units.SI.Time RVup(displayUnit="ms")=0.207
                                              " Duration of upstroke of RV activation";
  parameter Modelica.Units.SI.Time RVactdur0(displayUnit="ms")=0.529
                                                   " Duration of RV activation at HR 60";
    parameter Real RVdownpart=0.5
      " Duration of downstroke of RV activation as fraction of total duration";
    Real RVactdur=RVactdur0*((1/HRo)^(1/3))
      " Duration of RV activation, HR dependent, Fridericia";
    Real RVdown=RVactdur*RVdownpart " Duration of downstroke of RV activation";
  Real PI = Modelica.Constants.pi;

  // Cycle variables
    Modelica.Units.SI.Time t0(start = -0.1) "Beat time";
   Real CycleL=1/HRo;    // Cycle length
   Integer CycleCt(start = 0); // Cycle count,  = integer(floor((time)/CycleL))
   Real tcycle = time-t0; // Time within cycle
   Real tcycle2 = tcycle+CycleL;  // overlapping time into next cycle

  // Atria activation function
   Real Aact = if fixedActivation then fixedActValue else if (tcycle < Atact or tcycle > Atact+Aactdur) then 0 else
       (if (tcycle <= Atact+Aup) then (1-cos((tcycle-Atact)*PI/Aup))/2 else
       (if (tcycle <= Atact+Aactdur-Adown) then 1 else
        (1+cos((tcycle-(Atact+Aactdur-Adown))/(Adown)*PI))/2));

  //  Real LVact_n =        if (tcycle<LVtact) then
  //      (if  (tcycle2 > LVtact+LVactdur or CycleCt == 0) then 0 else
  //         noEvent(if (tcycle2 <= LVtact+LVup) then 1 else
  //           (if (tcycle2 <= LVtact+LVactdur-LVdown) then 2 else 3)))
  //  else
  //      (if  (tcycle > LVtact+LVactdur) then 4 else
  //        noEvent(if (tcycle <= LVtact+LVup) then 5 else
  //          (if (tcycle <= LVtact+LVactdur-LVdown) then 6 else 7)));

        // LV activation (including septum)
    Real LVact =    if fixedActivation then fixedActValue else  if (tcycle<LVtact) then
        (if  (tcycle2 > LVtact+LVactdur or CycleCt == 0) then 0 else
           (if (tcycle2 <= LVtact+LVup) then (sin((tcycle2-LVtact)*(PI/2)/LVup)) else
             (if (tcycle2 <= LVtact+LVactdur-LVdown) then 1 else
                (1+cos((tcycle2-(LVtact+LVactdur-LVdown))/(LVdown)*PI))/2)))
    else
        (if  (tcycle > LVtact+LVactdur) then 0 else
          (if (tcycle <= LVtact+LVup) then (sin((tcycle-LVtact)*(PI/2)/LVup)) else
            (if (tcycle <= LVtact+LVactdur-LVdown) then 1 else
                                                  (1+cos((tcycle-(LVtact+LVactdur-LVdown))/(LVdown)*PI))/2)));

  // RV activation
          Real RVact =        if fixedActivation then fixedActValue else if (tcycle<RVtact) then
     (if  (tcycle2 > RVtact+RVactdur or CycleCt == 0) then 0 else
       (if (tcycle2 <= RVtact+RVup) then (sin((tcycle2-RVtact)*(PI/2)/RVup)) else
         (if (tcycle2 <= RVtact+RVactdur-RVdown) then 1 else
                                                  (1+cos((tcycle2-(RVtact+RVactdur-RVdown))/(RVdown)*PI))/2)))
                          else
     (if  (tcycle > RVtact+RVactdur) then 0 else
       (if (tcycle <= RVtact+RVup) then (sin((tcycle-RVtact)*(PI/2)/RVup)) else
       (if (tcycle <= RVtact+RVactdur-RVdown) then 1 else
                                                  (1+cos((tcycle-(RVtact+RVactdur-RVdown))/(RVdown)*PI))/2)));

  // Pseudo-ECG
    Real pECG =  20*(-(LVact-0.5)^2+0.25) + 4*(-(Aact-0.5)^2+0.25);

  equation
    when time > pre(t0) + CycleL then
      CycleCt = pre(CycleCt) + 1;
      t0 = time;
    end when;
  //   HRo = if HR_stepUp then (60 + floor(time/10)*10)/60 else HRo0;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      experiment(
        StopTime=100,
        __Dymola_NumberOfIntervals=5000,
        __Dymola_Algorithm="Dassl"));
  end partialDriving_Olsen;

  model System_OlsenSIUnits "Base circulation for cotnrols"
    extends partialDriving_Olsen(t0(
          start=0));
    parameter Real rho=1060 " kg/m^3, // Density of blood";
    parameter Physiolibrary.Types.HydraulicCompliance C_ao(displayUnit="ml/mmHg")=
       8.1231668664085e-09   "mL/mmHg,  // Aortic compliance";
    parameter Physiolibrary.Types.HydraulicCompliance C_vc(displayUnit="ml/mmHg")=
       3.0002463033826e-07   "mL/mmHg,  // Vena cava compliance";
    parameter Physiolibrary.Types.HydraulicCompliance C_pa(displayUnit="ml/mmHg")=
       2.7752278306289e-08   "mL/mmHg,  // Pulmonary artery compliance";
    parameter Physiolibrary.Types.HydraulicCompliance C_pve(displayUnit="ml/mmHg")=
       1.1250923637685e-07   "mL/mmHg,  // Pulmonary veins compliance";
    //   Unstressed volumes
    parameter Physiolibrary.Types.Volume V0_ao=0
      "mL,   // Aorta unstressed volume";
    parameter Physiolibrary.Types.Volume V0_vc=0
      "mL,   // Vena cava unstressed volume";
    parameter Physiolibrary.Types.Volume V0_pa=0
      "mL,   // Pulmonary artery unstressed volume";
    parameter Physiolibrary.Types.Volume V0_pve=0
      "mL,   // Pulmonary veins unstressed volume";
    //   Inertances
    parameter Physiolibrary.Types.HydraulicInertance L_ao(displayUnit="mmHg.s2/ml")=
       1333223.87415   "mmHg*s^2*mL^(-1), // Aorta inertance";
    parameter Physiolibrary.Types.HydraulicInertance L_pa(displayUnit="mmHg.s2/ml")=
       1333223.87415   "mmHg*s^2*mL^(-1), // Pulmonary artery inertance";
    parameter Physiolibrary.Types.HydraulicInertance L_vc(displayUnit="mmHg.s2/ml")=
       1333223.87415   "mmHg*s^2*mL^(-1), // Vena cava inertance";
    parameter Physiolibrary.Types.HydraulicInertance L_pve(displayUnit="mmHg.s2/ml")=
       1333223.87415   "mmHg*s^2*mL^(-1), // Pulmonary veins inertance";

    //   Resistances
    parameter Physiolibrary.Types.HydraulicResistance R_vc(displayUnit="(mmHg.s)/ml")=
       133322.387415   "mmHg*s*mL^(-1), // Vena cava resistance";
    parameter Physiolibrary.Types.HydraulicResistance R_pve(displayUnit="(mmHg.s)/ml")=
       133322.387415   "mmHg*s*mL^(-1), // Pulmonary veins resistance";
    // Vascular resistances
    parameter Physiolibrary.Types.HydraulicResistance SVR(displayUnit="(mmHg.s)/ml")=
       160920121.60991   "mmHg*s*mL^(-1),  // Systemic vascular resistance";
    parameter Physiolibrary.Types.HydraulicResistance PVR(displayUnit="(mmHg.s)/ml")=
       10012511.294867   "mmHg*s*mL^(-1),  // Pulmonary vascular resistance";

    // Valves

    //   Valve areas
    parameter Physiolibrary.Types.Area AVAopen=0.00035
      "cm^2,  // Aortic valve area (when open)";
    parameter Physiolibrary.Types.Area AVAclosed(displayUnit="cm2")=1e-07
      "cm^2,  // (when closed)";
    parameter Physiolibrary.Types.Area MVAopen(displayUnit="cm2")=0.0004
      "cm^2,   // Mitral valve area (when open)";
    parameter Physiolibrary.Types.Area MVAclosed=1e-07 "cm^2,  // (when closed)";
    parameter Physiolibrary.Types.Area PVAopen=0.00035
      "cm^2,  // Pulmonary valve area (when open)";
    parameter Physiolibrary.Types.Area PVAclosed=1e-07 "cm^2,  // (when closed)";
    parameter Physiolibrary.Types.Area TVAopen=0.0004
      "cm^2,   // Tricuspid valve area (when open)";
    parameter Physiolibrary.Types.Area TVAclosed=1e-07 "cm^2,  // (when closed)";

    // Chambers
    //   Inertances
    parameter Physiolibrary.Types.HydraulicInertance L_lvot(displayUnit="mmHg.s2/ml")=
       133322.387415   "mmHg*s^2*mL^(-1),  // LVOT inertance";
    parameter Physiolibrary.Types.HydraulicInertance L_rvot(displayUnit="mmHg.s2/ml")=
       133322.387415   "mmHg*s^2*mL^(-1),  // RVOT inertance";
    parameter Physiolibrary.Types.HydraulicInertance L_la(displayUnit="mmHg.s2/ml")=
       13332.2387415   "mmHg*s^2*mL^(-1), // Left atrium inertance";
    parameter Physiolibrary.Types.HydraulicInertance L_ra(displayUnit="mmHg.s2/ml")=
       13332.2387415   "mmHg*s^2*mL^(-1), // Right atrium inertance";

    parameter Physiolibrary.Types.HydraulicInertance L_sept(displayUnit="mmHg.s2/ml")=
       33330.59685375   "mmHg*s^2*mL^(-1), // Septal displacement inertance";

    parameter Real Lm_ym=0.2 "g/m,   // Junction circle inertial component";
    parameter Real damp_ym=60  "N*s*m^(-1)*cm^(-1),   // Damping component";

    // LV - active myofiber properties (intercepts of assumed linear force-velocity relationship)
    parameter Modelica.Units.SI.Stress LV_smax(displayUnit="kPa")=134000
      "kPa, // Maximal active myofiber stress (at zero velocity)";
    parameter Real LV_dedtmax=-3
      "1/s, // Maximal myofiber strain rate (at zero afterload)";

    // LV - passive myofiber properties
    //   For extension above slack length and compression below slack length (assuming symmetric passive properties)
    parameter Physiolibrary.Types.Pressure LV_k(displayUnit="kPa")=4000
      "kPa,  // Parameters determining shape of compliance curve";
    parameter Real LV_a=11 "dimensionless,";
    //   Viscosity
    parameter Real LV_visc=1.5e3 "kPa*s, // Viscosity (stress per strain rate)";
    // Ref. Chung et al, J Mol Cell Cardiol. 2011 Sep; 51(3): 428–434.
    // They find titin-based viscosity approx 5.8 kPa at stretch rate of 1 s-1 (fig. 2)
    // in mice. A lower value gives more realistic results in this model.

    // RV - active myofiber properties
    parameter Physiolibrary.Types.Pressure RV_smax(displayUnit="kPa")=70000
      "kPa, // Maximal active myofiber stress";
    parameter Real RV_dedtmax=-3 "1/s, // Maximal myofiber strain rate";

    // RV - passive myofiber properties
    //   For extension above slack length and compression below slack length (assuming symmetric passive properties)
    parameter Physiolibrary.Types.Pressure RV_k(displayUnit="kPa")=4000   "kPa, ";
    parameter Real RV_a=11 "dimensionless,";
    //   Viscosity
    parameter Real RV_visc=1.5e3 "kPa*s, // viscosity (stress per strain rate)";

    //  Force-length relation
    parameter Real e0=-0.3
      "dimensionless,  // natural fiber strain associated with zero force";
    // force increases above this number (see formula in code)
    // LV chamber
    parameter Physiolibrary.Types.Volume LV_Vwall=0.0001015
      "cm^3, // Volume of LV wall (including septum)";
    parameter Physiolibrary.Types.Volume LV_V0=9.21e-05
      "cm^3,  // Unstressed volume of LV cavity";

    parameter Physiolibrary.Types.Volume LV_Vm0=LV_V0 + (LV_Vwall/2)
      ",  // Unstressed mid-wall volume of entire LV";
    Physiolibrary.Types.Volume r0_LV=(LV_Vm0*3/4/PI)^(1/3)
      ", // Unstressed mid-wall radius of spherical LV";

    // RV chamber
    parameter Physiolibrary.Types.Volume RV_Vwall=4.6e-05
      "cm^3, // Volume of RV wall (ref Katz et al., J Am Coll Cardiol. 1993 May;21(6):1475-81)";
    parameter Physiolibrary.Types.Volume RV_V0=0.000116
      "cm^3,  // Unstressed volume of RV cavity ";

    // Septum
    parameter Real Sept_ratio=0.333
      "dimensionless,  // Ratio of septal wall volume to LV wall volume";
    parameter Physiolibrary.Types.Volume Sept_Vwall=LV_Vwall*Sept_ratio
      "// Volume of septal wall";
    parameter Physiolibrary.Types.Volume LW_Vwall=LV_Vwall - Sept_Vwall
      "// Volume of left ventricle free wall";
    // (by definition, RV wall = RV free wall)

    // Spherical caps
    // For reference, unstressed values, a spherical shape of LV + septum is assumed
    Modelica.Units.SI.Height Sept_h0=2*Sept_ratio*r0_LV
      "  // Reference midwall height of spherical cap when ration btw cap area and spherical area is Sept_ratio";
    Real Sept_ym0=(2*Sept_h0*(r0_LV - 0.5*Sept_h0))^(1/2)
      " // Reference septum midwall junction radius";
    Physiolibrary.Types.Volume Sept_Vm0=(PI*Sept_h0^2)/3*(3*r0_LV - Sept_h0)
      "  // Reference septum midwall cap volume";
    Real Sept_Am0=2*PI*r0_LV*Sept_h0 "    // Reference septum midwall cap area";
    Real Sept_xm0;
    Real Sept_Cm0;

    Physiolibrary.Types.Volume LW_Vm0=-LV_V0 - 0.5*LW_Vwall - 0.5*Sept_Vwall +
        Sept_Vm0 " // Reference LV free wall midwall cap volume";
    Real LW_Am0=2*PI*r0_LV*(2*r0_LV - Sept_h0)
      "     // Reference LV free wall midwall cap area";
    Real LW_xm0;
    Real LW_Cm0;

    Physiolibrary.Types.Volume RW_Vm0=RV_V0 + 0.5*RV_Vwall + 0.5*Sept_Vwall +
        Sept_Vm0 "  // Reference RV free wall midwall cap volume";
    Real RW_h0 "  // Reference RV free wall midwall cap height";
    Real r0_RW ",   // Reference RV free wall midwall cap radius";

    Real RW_Am0 ",   // Reference RV free wall midwall cap area";
    Real RW_xm0;
    Real RW_Cm0;

    // Parameters, continued
    // LA
    parameter Physiolibrary.Types.Volume V0_la=2e-05
      "mL, // LA equilibrium volume";
    parameter Physiolibrary.Types.HydraulicCompliance C_la(displayUnit="ml/mmHg")=
       6.0004926067653e-08   "mL/mmHg, // LA compliance";
    parameter Physiolibrary.Types.Pressure P_maxla(displayUnit="mmHg")=799.93432449
      "mmHg, // Force of atrial contraction";

    // RA
    parameter Physiolibrary.Types.Volume V0_ra=2e-05
      "mL, // RA equilibrium volume";
    parameter Physiolibrary.Types.HydraulicCompliance C_ra(displayUnit="ml/mmHg")=
       6.0004926067653e-08   "mL/mmHg, // RA compliance";
    parameter Physiolibrary.Types.Pressure P_maxra(displayUnit="mmHg")=533.28954966
      "mmHg, // Force of RA contraction";

    // Pericardium
    parameter Physiolibrary.Types.Volume V0_peri=0.000599
      "mL,  // Pericardium equilibrium volume";
    parameter Physiolibrary.Types.Pressure k_peri(displayUnit="mmHg")=133.322387415
      "mmHg,  // Compliance curve constant";
    parameter Real a_peri=0.017e6 "mL^(-1), // Compliance curve constant";
    parameter Physiolibrary.Types.Volume V_perifl=1e-05
      "mL,  // Pericardial fluid volume";

    // Initial conditions
    //
    parameter Physiolibrary.Types.Volume V_lvinit=LV_V0
      ",                        // LV initial volume";
    parameter Physiolibrary.Types.Volume V_rvinit=RV_V0
      ",                        // RV initial volume";
    parameter Physiolibrary.Types.Volume V_aoinit=0.0001
      "mL ,                        // Aorta initial volume";
    parameter Physiolibrary.Types.Volume V_vcinit=0.000123
      "mL ,                        // Vena cava initial volume";
    parameter Physiolibrary.Types.Volume V_painit=7e-05
      "mL ,                        // PA initial volume";
    parameter Physiolibrary.Types.Volume V_pveinit=0.00015
      "mL ,                        // PVe initial volume";
    parameter Physiolibrary.Types.Volume V_rainit=5e-05
      "mL ,                        // RA initial volume";
    parameter Physiolibrary.Types.Volume V_lainit=8e-05
      "mL,                        // LA initial volume";

    // -----------------------------------------------------------------------------
    //   VARIABLES
    // -----------------------------------------------------------------------------
    //  Pressures
    Physiolibrary.Types.Pressure P_lvw " mmHg,   // LV pressure at wall";

    Physiolibrary.Types.Pressure P_rvw " mmHg,   // RV pressure at wall";

    Physiolibrary.Types.Pressure P_la "  mmHg,   // LA pressure";
    Physiolibrary.Types.Pressure P_ra "  mmHg,   // RA pressure";

    Physiolibrary.Types.Pressure P_ao "  mmHg,   // Aortic pressure";
    Physiolibrary.Types.Pressure P_vc "  mmHg,   // Vena cava pressure";
    Physiolibrary.Types.Pressure P_pa "  mmHg,   // Pulmonary arterial pressure";
    Physiolibrary.Types.Pressure P_pve " mmHg,   // Pulmonary veins pressure";

    Physiolibrary.Types.Pressure P_peri " mmHg,   // Pericardial pressure";

    //  Volumes
    Physiolibrary.Types.Volume V_lv "  mL,   // LV volume";
    Physiolibrary.Types.Volume V_rv "  mL,   // RV volume";
    Physiolibrary.Types.Volume V_la "  mL,   // LA volume";
    Physiolibrary.Types.Volume V_ra "  mL,   // RA volume";
    Physiolibrary.Types.Volume V_ao "  mL,   // Aortic volume";
    Physiolibrary.Types.Volume V_vc "  mL,   // Vena cava volume";
    Physiolibrary.Types.Volume V_pa "  mL,   // Pulmonary artery volume";
    Physiolibrary.Types.Volume V_pve " mL,   // Pulmonary veins volume";

    Physiolibrary.Types.Volume V_peri " mL,   // Total pericardial volume";

    //  Flows
    Physiolibrary.Types.VolumeFlowRate Q_vc
      "  mL/s,   // Flow from vena cava to RA";
    Physiolibrary.Types.VolumeFlowRate Q_tv
      "  mL/s,   // Flow through tricuspid valve";
    Physiolibrary.Types.VolumeFlowRate Q_pva
      " mL/s,   // Flow through pulmonary valve";
    Physiolibrary.Types.VolumeFlowRate Q_pa
      "  mL/s,   // Flow from PA to pulmonary veins";
    Physiolibrary.Types.VolumeFlowRate Q_pve
      " mL/s,   // Flow from pulmonary veins to LA";
    Physiolibrary.Types.VolumeFlowRate Q_mv
      "  mL/s,   // Flow through mitral valve";
    Physiolibrary.Types.VolumeFlowRate Q_av
      "  mL/s,   // Flow through aortic valve";
    Physiolibrary.Types.VolumeFlowRate Q_ao
      "  mL/s,   // Flow from aorta to vena cava";

    // Valve areas
    Physiolibrary.Types.Area AVA "  cm^2,   // Time-dependent aortic valve area";
    Physiolibrary.Types.Area MVA "  cm^2,   // Time-dependent mitral valve area";
    Physiolibrary.Types.Area PVA
      "  cm^2,   // Time-dependent pulmonary valve area";
    Physiolibrary.Types.Area TVA
      "  cm^2,   // Time-dependent tricuspid valve area";

    //  Spherical cap variables
    Modelica.Units.SI.Radius Sept_ym "  cm,   // Septum midwall junction radius";

    Physiolibrary.Types.Length Sept_xm
      "  cm,   // Septum midwall axial distance from origin";
    Physiolibrary.Types.Volume Sept_Vm "  cm^3,   // Septum midwall cap volume";
    Physiolibrary.Types.Area Sept_Am "  cm^2,   // Septum midwall cap area";
    Real Sept_Cm "  cm^-1,    // Septum midwall cap curvature (1/r)";
    Modelica.Units.SI.SurfaceTension Sept_Tm "  N/m,   // Septum midwall tension";
    Modelica.Units.SI.SurfaceTension Sept_Tx
      "  N/m,   // Septum x-component tension at junction (axial)";
    Modelica.Units.SI.SurfaceTension Sept_Ty
      "  N/m,   // Septum y-component tension at junction (radial)";

    Physiolibrary.Types.Pressure P_Septtrans
      "  kPa,   // Septum transmural pressure (signed)";

    Real Sept_e "  dimensionless,  // Myofiber natural strain";
    Modelica.Units.SI.Frequency Sept_dedt "  s^(-1),   // Myofiber strain rate";
    Physiolibrary.Types.Pressure Sept_s_act "  kPa,   // Myofiber active stress";
    Physiolibrary.Types.Pressure Sept_s_ela
      "  kPa,   // Myofiber passive stress from elasticity";
    Physiolibrary.Types.Pressure Sept_s_vis
      "  kPa,   // Myofiber passive stress from viscous resistance to motion";
    Physiolibrary.Types.Pressure Sept_s "  kPa,   // Total myofiber stress";
    Real Sept_z "  dimensionless,  // Dimensionless curvature parameter";

    Physiolibrary.Types.Length LW_xm
      "  cm,   // LV free wall midwall axial distance from origin";
    Physiolibrary.Types.Volume LW_Vm
      "  cm^3,   // LV free wall midwall cap volume";
    Physiolibrary.Types.Area LW_Am "  cm^2,   // LV free wall midwall cap area";
    Physiolibrary.Types.Area LW_Cm
      "  cm^-1,    // LV free wall midwall cap curvature (1/r)";
    Modelica.Units.SI.SurfaceTension LW_Tm "  N/m,   // LV midwall tension";
    Modelica.Units.SI.SurfaceTension LW_Tx
      "  N/m,   // LV x-component tension at junction (axial)";
    Modelica.Units.SI.SurfaceTension LW_Ty
      "  N/m,   // LV y-component tension at junction (radial)";

    Physiolibrary.Types.Pressure P_LWtrans
      "  kPa,   // LV free wall transmural pressure (signed)";

    Real LW_e "   dimensionless,  // Myofiber natural strain";
    Modelica.Units.SI.Frequency LW_dedt "  s^(-1),   // Myofiber strain rate";
    Physiolibrary.Types.Pressure LW_s_act "  kPa,   // Myofiber active stress";
    Physiolibrary.Types.Pressure LW_s_ela
      "  kPa,   // Myofiber passive stress from elasticity";
    Physiolibrary.Types.Pressure LW_s_vis
      "  kPa,   // Myofiber passive stress from viscous resistance to motion";
    Physiolibrary.Types.Pressure LW_s "   kPa,   // Total myofiber stress";
    Real LW_z "   dimensionless,  // Dimensionless curvature parameter";

    Physiolibrary.Types.Length RW_xm
      "  cm,   // RV free wall midwall axial distance from origin";
    Physiolibrary.Types.Volume RW_Vm
      "  cm^3,   // RV free wall midwall cap volume";
    Physiolibrary.Types.Area RW_Am "  cm^2,   // RV free wall midwall cap area";
    Real RW_Cm "  cm^-1,    // RV free wall midwall cap curvature (1/r)";
    Modelica.Units.SI.SurfaceTension RW_Tm "  N/m,   // RV midwall tension";
    Modelica.Units.SI.SurfaceTension RW_Tx
      "  N/m,   // RV x-component tension at junction (axial)";
    Modelica.Units.SI.SurfaceTension RW_Ty
      "  N/m,   // RV y-component tension at junction (radial)";

    Real P_RWtrans "  kPa,   // RV free wall transmural pressure (signed)";

    Real RW_e "   dimensionless,  // Myofiber natural strain";
    Modelica.Units.SI.Frequency RW_dedt "  s^(-1),   // Myofiber strain rate";
    Physiolibrary.Types.Pressure RW_s_act "  kPa,   // Myofiber active stress";
    Physiolibrary.Types.Pressure RW_s_ela
      "  kPa,   // Myofiber passive stress from elasticity";
    Physiolibrary.Types.Pressure RW_s_vis
      "  kPa,   // Myofiber passive stress from viscous resistance to motion";
    Physiolibrary.Types.Pressure RW_s "   kPa,   // Total myofiber stress";
    Real RW_z "   dimensionless,  // Dimensionless curvature parameter";

    Modelica.Units.SI.SurfaceTension Tx_total
      "  N/m,   // Sum tension at junction circle in x-direction (axial)";
    Modelica.Units.SI.SurfaceTension Ty_total
      "  N/m,   // Sum tension at junction circle in y-direction (radial)";

    Physiolibrary.Types.VolumeFlowRate Q_sept
      "  mL/s,   // Flow into septal displacement volume";

    Physiolibrary.Types.Velocity d_Sept_ym
      "  cm/s;   // Change in junction radius diameter";

  initial equation
    // -----------------------------------------------------------------------------
    //   INITIAL CONDITIONS
    // -----------------------------------------------------------------------------

    //        Initial Conditions
    //
    //           when(t=t.min) {
    V_ao = V_aoinit;
    V_vc = V_vcinit;
    V_pa = V_painit;
    V_pve = V_pveinit;
    V_lv = V_lvinit;
    V_rv = V_rvinit;
    V_la = V_lainit;
    V_ra = V_rainit;

    Q_tv = 0;
    Q_pva = 0;
    Q_pa = 0;
    Q_vc = 0;
    Q_pve = 0;
    Q_mv = 0;
    Q_av = 0;
    Q_ao = 0;

    Q_sept = 0;

    d_Sept_ym = 0;

    Sept_Vm = Sept_Vm0;

    Sept_ym = Sept_ym0;

    // State variables
    //         AVclosed        = 1;
    //         MVclosed        = 0;
    //         PVclosed        = 1;
    //         TVclosed        = 0;
    //
    //           }

  equation
    HRo = if HR_stepUp then (60 + floor(time/stepUpLength) *10)/60 else HRo0;
    // Implicit equations, spherical caps
    r0_RW = (Sept_ym0^2 + RW_h0^2)/(2*RW_h0);
    RW_Vm0 = (PI*RW_h0^2)/3*(3*r0_RW - RW_h0);
    RW_Am0 = 2*PI*r0_RW*RW_h0;

    LW_Vm0 = PI/6*LW_xm0*(LW_xm0^2 + 3*Sept_ym0^2);
    // Yields LW_xm0
    LW_Cm0 = 2*LW_xm0/(LW_xm0^2 + Sept_ym0^2);

    Sept_Vm0 = PI/6*Sept_xm0*(Sept_xm0^2 + 3*Sept_ym0^2);
    Sept_Cm0 = 2*Sept_xm0/(Sept_xm0^2 + Sept_ym0^2);

    RW_Vm0 = PI/6*RW_xm0*(RW_xm0^2 + 3*Sept_ym0^2);
    RW_Cm0 = 2*RW_xm0/(RW_xm0^2 + Sept_ym0^2);

    // -----------------------------------------------------------------------------
    //   SYSTEM OF EQUATIONS
    // -----------------------------------------------------------------------------

    // Ordinary differential equations

    // Valve flows
    // Tricuspid valve
    der(Q_tv) = (P_ra - P_rvw - ((Q_tv*abs(Q_tv)/(TVA^2))*rho/2))/L_ra;
    // Pulmonary valve
    der(Q_pva) = (P_rvw - P_pa - ((Q_pva*abs(Q_pva)/(PVA^2))*rho/2))/L_rvot;
    // Mitral valve
    der(Q_mv) = (P_la - P_lvw - ((Q_mv*abs(Q_mv)/(MVA^2))*rho/2))/L_la;
    // Aortic valve
    der(Q_av) = (P_lvw - P_ao - ((Q_av*abs(Q_av)/(AVA^2))*rho/2))/L_lvot;

    // Vessel flows
    // Vena cava
    der(Q_vc) = (P_vc - P_ra - Q_vc*R_vc)/L_vc;
    // Pulmonary artery
    der(Q_pa) = (P_pa - P_pve - Q_pa*PVR)/L_pa;
    // Pulmonary veins
    der(Q_pve) = (P_pve - P_la - Q_pve*R_pve)/L_pve;
    // Aorta
    der(Q_ao) = (P_ao - P_vc - Q_ao*SVR)/L_ao;

    // Conservation of mass equations

    // der(V_vc)  = Q_ao - Q_vc;
    // der(V_ra)  = 0;//!Q_vc - Q_tv;
    // der(V_rv)  = 0;//!Q_tv - Q_pva;
    // der(V_pa)  = Q_pva - Q_pa;
    // der(V_pve)  = Q_pa - Q_pve;
    // der(V_la)  = 0;//!Q_pve - Q_mv;
    // der(V_lv)  =  0;//!Q_mv - Q_av;
    // der(V_ao)  = Q_av - Q_ao;
    der(V_vc) = Q_ao - Q_vc;
    der(V_ra) = Q_vc - Q_tv;
    der(V_rv) = Q_tv - Q_pva;
    der(V_pa) = Q_pva - Q_pa;
    der(V_pve) = Q_pa - Q_pve;
    der(V_la) = Q_pve - Q_mv;
    der(V_lv) = Q_mv - Q_av;
    der(V_ao) = Q_av - Q_ao;

    // Change in valve area with state
    AVA = if (Q_av < 0) then AVAclosed else AVAopen;
    MVA = if (Q_mv < 0) then MVAclosed else MVAopen;
    PVA = if (Q_pva < 0) then PVAclosed else PVAopen;
    TVA = if (Q_tv < 0) then TVAclosed else TVAopen;

    // Spherical caps geometries
    // Displacement flow - interventricular septum
    der(Q_sept) = -(P_LWtrans + P_Septtrans + P_RWtrans)/L_sept;
    // Flow acceleration into septal spherical cap

    der(Sept_Vm) = Q_sept;
    // Change in septal midwall volume (positive when curved towards RV)

    LW_Vm = -V_lv - 0.5*LW_Vwall - 0.5*Sept_Vwall + Sept_Vm;
    // LV free wall spherical cap volume - adjusted according to septal displacement
    RW_Vm = V_rv + 0.5*RV_Vwall + 0.5*Sept_Vwall + Sept_Vm;
    // RV free wall spherical cap volume - same

    // Geometry changes of septal junction circle

    der(d_Sept_ym) = -(Ty_total + d_Sept_ym*damp_ym)/Lm_ym;
    // Acceleration in junction radius caused by summed tension forces + damping
    der(Sept_ym) = d_Sept_ym;
    // Change in junction radius

    // Solve for xm - cap midwall height (positive in direction of RV)
    LW_Vm = PI/6*LW_xm*(LW_xm^2 + 3*Sept_ym^2);
    // Given Vm and ym, xm can be found
    RW_Vm = PI/6*RW_xm*(RW_xm^2 + 3*Sept_ym^2);
    Sept_Vm = PI/6*Sept_xm*(Sept_xm^2 + 3*Sept_ym^2);

    // Calculate Am - cap midwall area
    LW_Am = PI*(LW_xm^2 + Sept_ym^2);
    RW_Am = PI*(RW_xm^2 + Sept_ym^2);
    Sept_Am = PI*(Sept_xm^2 + Sept_ym^2);

    // Calculate Cm - cap midwall curvature (1/r)
    LW_Cm = 2*LW_xm/(LW_xm^2 + Sept_ym^2);
    RW_Cm = 2*RW_xm/(RW_xm^2 + Sept_ym^2);
    Sept_Cm = 2*Sept_xm/(Sept_xm^2 + Sept_ym^2);

    // Tx balance
    Tx_total = LW_Tx + Sept_Tx + RW_Tx;
    // Ty balance
    Ty_total = LW_Ty + Sept_Ty + RW_Ty;

    // Fiber strain
    // Relation btw fiber strain and cavity + wall volume
    // From Lumens et al., 2009, Annals of Biomedical Engineering.
    LW_z = (3*LW_Cm*LW_Vwall)/(2*LW_Am);
    // Curvature ratio
    LW_e = 0.5*log(LW_Am/LW_Am0) - (1/12)*LW_z^2 - 0.019*LW_z^4;
    // Natural strain. Approximation to solution for curved wall segment (Lumens 2009)
    LW_dedt = der(LW_e);
    // Strain rate

    RW_z = (3*RW_Cm*RV_Vwall)/(2*RW_Am);
    // Curvature ratio
    RW_e = 0.5*log(RW_Am/RW_Am0) - (1/12)*RW_z^2 - 0.019*RW_z^4;
    // Natural strain. Approximation to solution for curved wall segment (Lumens 2009)
    RW_dedt = der(RW_e);
    // Strain rate

    Sept_z = (3*Sept_Cm*Sept_Vwall)/(2*Sept_Am);
    // Curvature ratio
    Sept_e = 0.5*log(Sept_Am/Sept_Am0) - (1/12)*Sept_z^2 - 0.019*Sept_z^4;
    // Natural strain. Approximation to solution for curved wall segment (Lumens 2009)
    Sept_dedt = der(Sept_e);
    // Strain rate

    // Fiber stress, active component (applying force-velocity and force-length relation)
    // A linear relation btw ratio of developed force vs. max force
    // and ratio of strain rate vs. max strain rate is assumed.
    // This seems reasonable for the relationship at physiological loads
    // (force 15-100% of peak force). A hyperbolic relationship
    // breaks down above 70% of peak force.
    // Note that max strain rate for this assumption is different from
    // and lower than Vmax estimated when using the hyperbolic Hill equation.
    // ref. Circ Res. 1991 Feb;68(2):588-96. Sarcomere dynamics in cat cardiac trabeculae. de Tombe PP1, ter Keurs HE.
    LW_s_act = LVact*(if (LW_dedt > 0) then (LV_smax*tanh((LW_e - e0)*3)) else (
      if (LW_dedt < LV_dedtmax) then 0 else (LV_smax - LV_smax*(LW_dedt/
      LV_dedtmax))*tanh((LW_e - e0)*3)));

    Sept_s_act = LVact*(if (Sept_dedt >= 0) then (LV_smax*tanh((Sept_e - e0)*3))
       else (if (Sept_dedt < LV_dedtmax) then 0 else (LV_smax - LV_smax*(
      Sept_dedt/LV_dedtmax))*tanh((Sept_e - e0)*3)));

    RW_s_act = RVact*(if (RW_dedt >= 0) then (RV_smax*tanh((RW_e - e0)*3)) else (
      if (RW_dedt < RV_dedtmax) then 0 else (RV_smax - RV_smax*(RW_dedt/
      RV_dedtmax))*tanh((RW_e - e0)*3)));

    // Fiber stress, passive components
    // Elastic stress modeled to be symmetric on both sides of slack length
    // Viscous force modeled to be proportional to strain rate

    LW_s_ela = if (LW_e >= 0) then (LV_k*(exp(LV_a*LW_e) - 1)) else -(LV_k*(exp(
      LV_a*(-LW_e)) - 1));

    LW_s_vis = LW_dedt*LV_visc;

    Sept_s_ela = if (Sept_e >= 0) then (LV_k*(exp(LV_a*Sept_e) - 1)) else -(LV_k*(
      exp(LV_a*(-Sept_e)) - 1));

    Sept_s_vis = Sept_dedt*LV_visc;

    RW_s_ela = if (RW_e >= 0) then (RV_k*(exp(RV_a*RW_e) - 1)) else -(RV_k*(exp(
      RV_a*(-RW_e)) - 1));

    RW_s_vis = RW_dedt*RV_visc;

    // Fiber stress, total (sum of active and passive components)

    LW_s = LW_s_act + LW_s_ela + LW_s_vis;
    Sept_s = Sept_s_act + Sept_s_ela + Sept_s_vis;
    RW_s = RW_s_act + RW_s_ela + RW_s_vis;

    // Wall tension (Approximation, Lumens 2009)
    LW_Tm = (LW_Vwall*LW_s)/(2*LW_Am)*(1 + LW_z^2/3 + LW_z^4/5);
    LW_Tx = LW_Tm*(2*LW_xm*Sept_ym)/(LW_xm^2 + Sept_ym^2);
    LW_Ty = LW_Tm*(-(LW_xm^2) + Sept_ym^2)/(LW_xm^2 + Sept_ym^2);

    Sept_Tm = (Sept_Vwall*Sept_s)/(2*Sept_Am)*(1 + Sept_z^2/3 + Sept_z^4/5);
    Sept_Tx = Sept_Tm*(2*Sept_xm*Sept_ym)/(Sept_xm^2 + Sept_ym^2);
    Sept_Ty = Sept_Tm*(-(Sept_xm^2) + Sept_ym^2)/(Sept_xm^2 + Sept_ym^2);

    RW_Tm = (RV_Vwall*RW_s)/(2*RW_Am)*(1 + RW_z^2/3 + RW_z^4/5);
    RW_Tx = RW_Tm*(2*RW_xm*Sept_ym)/(RW_xm^2 + Sept_ym^2);
    RW_Ty = RW_Tm*(-(RW_xm^2) + Sept_ym^2)/(RW_xm^2 + Sept_ym^2);

    // Axial transmural pressure (in x-direction)
    P_LWtrans = 2*LW_Tx/Sept_ym;
    P_Septtrans = 2*Sept_Tx/Sept_ym;
    P_RWtrans = 2*RW_Tx/Sept_ym;

    // LV pressure
    P_lvw = -P_LWtrans + P_peri;

    // RV pressure

    P_rvw = P_RWtrans + P_peri;

    // Atrial pressures
    // Passive component + active component
    // Linear compliance curve
    P_la = (V_la - V0_la)/C_la + Aact*P_maxla + P_peri;
    P_ra = (V_ra - V0_ra)/C_ra + Aact*P_maxra + P_peri;

    // Vessels pressures

    P_ao = (V_ao - V0_ao)/C_ao;
    P_vc = (V_vc - V0_vc)/C_vc;
    P_pa = (V_pa - V0_pa)/C_pa;
    P_pve = (V_pve - V0_pve)/C_pve;

    // Pericardial volume and pressure
    V_peri = V_lv + V_rv + V_ra + V_la + LV_Vwall + RV_Vwall + V_perifl;
    P_peri = k_peri*(exp((V_peri - V0_peri)*a_peri) - 1);

    annotation (experiment(
        StopTime=20,
        Interval=0.001,
        Tolerance=1e-07,
        __Dymola_Algorithm="Dassl"));
  end System_OlsenSIUnits;

  model System_OlsenSIUnits_PAH
    "Basic pulmonary arterial hypertension parametrization"
    extends System_OlsenSIUnits(
      HRo0=1.1333333333333,
      V_lainit=0.00011,
      V_rainit=7e-05,
      V_pveinit=0.00021,
      V_vcinit=0.000315,
      V0_peri=0.000611,
      RV_V0=0.0001066,
      RV_Vwall=6.77e-05,
      LV_V0=7.05e-05,
      LV_Vwall=9.37e-05,
      RV_smax=127000,
      LV_smax=107000,
      PVR=80966685.877129,
      SVR=167186273.81841,
      C_pa=7.265096423641e-09,
      C_ao=6.9275687145105e-09,
      RVactdur0=0.585,
      RVdelay0=0.186,
      LVactdur0=0.564,
      LVup=0.255,
      LVdelay0=0.17);
  end System_OlsenSIUnits_PAH;

  model System_OlsenSIUnits_PAHEx
    "Parametrization to pulmonary arterial hypertension with exercise"
    extends System_OlsenSIUnits(
      HRo0(displayUnit="1/min") = 1.7833333333333,
      V_lainit=0.00011,
      V_rainit=0.0001,
      V_pveinit=0.00025,
      V_painit=0.0001,
      V_vcinit=0.000355,
      V0_peri=0.000631,
      P_maxra=799.93432449,
      P_maxla=1199.901486735,
      RV_V0=0.000116,
      RV_Vwall=7.3e-05,
      LV_V0=7.1e-05,
      LV_Vwall=9.95e-05,
      RV_smax=202000,
      LV_smax=161000,
      L_ra=13332.2387415,
      L_la=13332.2387415,
      LVdelay0=0.196,
      LVup=0.157,
      LVactdur0=0.477,
      RVdelay0=0.185,
      RVup=0.123,
      RVactdur0=0.53,
      C_ao(displayUnit="ml/mmHg") = 7.0490786897975e-09,
      C_vc(displayUnit="ml/mmHg") = 1.5001231516913e-07,
      C_pa(displayUnit="ml/mmHg") = 7.0490786897975e-09,
      SVR(displayUnit="(mmHg.s)/ml") = 115110549.29411,
      PVR(displayUnit="(mmHg.s)/ml") = 86219587.941281);
  end System_OlsenSIUnits_PAHEx;
  annotation (uses(
      Physiolibrary(version="2.4.1"),
      Modelica(version="4.0.0"),
      ADAN_main(version="1.1")));
end OlsenCirculation;
