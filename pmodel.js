
//----------------------------------------------------------------------- 
// Example values
//----------------------------------------------------------------------- 
var co2   = ee.Image(376.0)
//var elv   = ee.Image(450.0)
var elv   = ee.Image('CGIAR/SRTM90_V4')
var mtemp = [0.4879904, 6.1999985, 7.4999870, 9.6999003, 13.1999913, 19.6999227, 18.6000030, 18.0999577, 13.8999807, 10.7000307, 7.2999217, 4.4999644]
var mvpd  = [113.0432, 338.4469, 327.1185, 313.8799, 247.9747, 925.9489, 633.8551, 497.6772, 168.7784, 227.1889, 213.0142, 172.6035]
var mppfd = [223.8286, 315.2295, 547.4822, 807.4035, 945.9020, 1194.1227, 1040.5228, 1058.4161, 814.2580, 408.5199, 268.9183, 191.4482]
var fapar = ee.Image(1.0)

//-----------------------------------------------------------------------
// Define PFT-specific parameters. Here for C3 plants.
//-----------------------------------------------------------------------
var params_pft_gpp = {
  kphio       : 0.0744,
  beta        : 146.0,
  rd_to_vcmax : 0.015
}

//-----------------------------------------------------------------------
// PFT-independent parameters, therefore hard-wired.
//-----------------------------------------------------------------------
var kPo = 101325.0        // standard atmosphere, Pa (Allen, 1973)
var kTo = 25.0            // base temperature, deg C (Prentice, unpublished)
var temp0 = 0.0           // temperature below which all quantities are zero (deg C)
var temp1 = 12.0          // temperature above which ramp function is 1.0 (linear between temp0 and temp1) (deg C)

// Metabolic N ratio (N per unit Vcmax)
// Reference: Harrison et al., 2009, Plant, Cell and Environment; Eq. 3
var mol_weight_rubisco    = 5.5e5    // molecular weight of Rubisco, (g R)(mol R)-1
var n_conc_rubisco        = 1.14e-2  // N concentration in rubisco, (mol N)(g R)-1
var cat_turnover_per_site = 2.33     // catalytic turnover rate per site at 25 deg C, (mol CO2)(mol R sites)-1; use 2.33 instead of (3.5) as not all Rubisco is active (see Harrison et al., 2009)  
var cat_sites_per_mol_R   = 8.0      // number of catalytic sites per mol R, (mol R sites)(mol R)-1

// Metabolic N ratio (= 336.3734 mol N s (mol CO2)-1 )
var n_v = mol_weight_rubisco * n_conc_rubisco / ( cat_turnover_per_site * cat_sites_per_mol_R )


//-----------------------------------------------------------------------
// Execute P-model with specified arguments
//-----------------------------------------------------------------------
var result = pmodel( mtemp[0], mvpd[0], co2, elv, mppfd[0], fapar, "C3_full" )
print(result)

var out_keys = result.keys();
var out_model = ee.Image(0);

  var i_key = out_keys.get(1)
  var out_model = ee.Image(result.get(i_key));

for (var ii = 1; ii < 19; ii++) { 
  var i_key = out_keys.get(ii)
  var imTemp = ee.Image(result.get(i_key));
  var out_model = out_model.addBands(imTemp)
 
}

var out_model = out_model.rename(out_keys)

Map.addLayer(out_model.select('gpp'), {min:5,max:12}, 'gpp')


//-----------------------------------------------------------------------
// P-model
//-----------------------------------------------------------------------
function pmodel( tc, vpd, co2, elv, ppfd, fapar, method ) {
  
  var tc = ee.Image(tc);
  
  // absorbed photosynthetically active radiation (mol/m2)
  var ppfdabs = fapar.multiply(ee.Image(ppfd));
  
  // atmospheric pressure as a function of elevation (Pa)
  var patm = calc_patm( elv );

  // ambient CO2 partial pression (Pa)
  var ca = co2_to_ca( co2, patm );
  
  // photorespiratory compensation point - Gamma-star (Pa)
  var gstar = calc_gstar( tc );
  // print(gstar,'gstar')
   
  // Michaelis-Menten coef. (Pa)
  var kmm  = calc_k( tc, patm );
  // print(kmm,'kmm')

  // viscosity correction factor = viscosity( temp, press )/viscosity( 25 degC, 1013.25 Pa) 
  var ns      = calc_viscosity_h2o( tc, patm );  // Pa s 
  var ns25    = calc_viscosity_h2o( ee.Image(kTo), ee.Image(kPo) );  // Pa s 
  var ns_star = ns.divide(ns25);                       // (unitless)
 
  
  if (method=="approx"){
    //-----------------------------------------------------------------------
    // A. APPROXIMATIVE METHOD
    //-----------------------------------------------------------------------
    var out_lue = lue_approx( tc, vpd, elv, ca, gstar, ns, kmm );
              
  } else if (method=="C3_simpl"){
    //-----------------------------------------------------------------------
    // B.1 SIMPLIFIED FORMULATION 
    //-----------------------------------------------------------------------
    var out_lue = lue_vpd_c3_simpl( kmm, gstar, ns, ca, vpd );

 } else if (method=="C3_full"){
    //-----------------------------------------------------------------------
    // B.2 FULL FORMULATION
    //-----------------------------------------------------------------------
    var out_lue = lue_vpd_c3_full( kmm, gstar, ns_star, ca, vpd );

  } else if (method=="C4"){
    //-----------------------------------------------------------------------
    // B.2 FULL FORMULATION
    //-----------------------------------------------------------------------
    var out_lue = lue_c4();

  }

  
  // LUE-functions return m, n, and chi
  var chi = out_lue.chi;

  // print('m', out_lue.m )
  // print('chi', out_lue.chi )

  //-----------------------------------------------------------------------
  // Calculate function return variables
  //-----------------------------------------------------------------------
  // GPP per unit ground area is the product of the intrinsic quantum 
  // efficiency, the absorbed PAR, the function of alpha (drought-reduction),
  // and 'm'
  var mprime = calc_mprime( out_lue.m );
  
  // Light use efficiency (assimilation rate per unit absorbed light)
  var lue = ee.Image(params_pft_gpp.kphio).multiply(mprime) ; // in mol CO2 m-2 s-1 / (mol light m-2 s-1)
  

  // Gross primary productivity = ecosystem-level assimilation rate (per unit ground area)
  var assim = ppfdabs.multiply(lue); // in mol CO2 m-2 s-1
  
  // Leaf-level assimilation rate (per unit leaf area), representative for top-canopy leaves
  var assim_unitfapar = ee.Image(ppfd).multiply(lue);  // in mol m-2 s-1

  // leaf-internal CO2 partial pressure (Pa)
  var ci = chi.multiply(ca);
 
  

  // stomatal conductance to H2O, expressed per unit absorbed light
  var gs_unitiabs = ee.Image(1.6).multiply(lue).multiply(patm).divide( ca.subtract(ci) );


  // Vcmax per unit ground area is the product of the intrinsic quantum 
  // efficiency, the absorbed PAR, and 'n'
  var vcmax = ppfdabs.multiply(ee.Image(params_pft_gpp.kphio)).multiply(ee.Image(out_lue.n));
  // print( 'ppfdabs', ppfdabs )
  // print( 'params_pft_gpp.kphio', params_pft_gpp.kphio )
  // print( 'out_lue.n', out_lue.n )
  // print( 'vcmax', vcmax )

  // Vcmax normalised per unit fAPAR (assuming fAPAR=1)
  var vcmax_unitfapar = ee.Image(ppfd).multiply(ee.Image(params_pft_gpp.kphio)).multiply(ee.Image(out_lue.n));

  // Vcmax normalised per unit absorbed PPFD (assuming ppfdabs=1)
  var vcmax_unitiabs = ee.Image(params_pft_gpp.kphio).multiply(ee.Image(out_lue.n)); 

  // Vcmax25 (vcmax normalized to 25 deg C)
  var factor25_vcmax    = calc_vcmax25( ee.Image(1.0), tc );
  
  var vcmax25           = factor25_vcmax.multiply(vcmax);
  var vcmax25_unitfapar = factor25_vcmax.multiply(vcmax_unitfapar);
  var vcmax25_unitiabs  = factor25_vcmax.multiply(vcmax_unitiabs);
  // print( 'vcmax25', vcmax25 )

  // Dark respiration
  var rd = ee.Image(params_pft_gpp.rd_to_vcmax).multiply(vcmax);
 

  
  // Dark respiration per unit fAPAR (assuming fAPAR=1)
  var rd_unitfapar = ee.Image(params_pft_gpp.rd_to_vcmax).multiply(vcmax_unitfapar);

  // Dark respiration per unit absorbed PPFD (assuming ppfdabs=1)
  var rd_unitiabs = ee.Image(params_pft_gpp.rd_to_vcmax).multiply(vcmax_unitiabs);

  // active metabolic leaf N (canopy-level), mol N/m2-ground (same equations as for nitrogen content per unit leaf area, gN/m2-leaf)
  var actnv = vcmax25.multiply(ee.Image(n_v));
  var actnv_unitfapar = vcmax25_unitfapar.multiply(ee.Image(n_v));
  var actnv_unitiabs  = vcmax25_unitiabs.multiply(ee.Image(n_v));
  

  // Construct derived type for output
  var out_pmodel = ee.Dictionary({
    gpp:             assim,
    gstar:           gstar,
    chi:             chi,
    ci:              co2.multiply(chi),  // return value 'ci:' is used for output in units of ppm. ,
    ca:              ca,
    iwue:            ( ca.subtract(ci) ).divide( ee.Image(1.6).multiply(patm) ),
    gs_unitiabs:     gs_unitiabs,
    vcmax:           vcmax,
    vcmax25:         vcmax25,
    vcmax_unitfapar: vcmax_unitfapar,
    vcmax_unitiabs:  vcmax_unitiabs,
    factor25_vcmax:  factor25_vcmax,
    rd:              rd,
    rd_unitfapar:    rd_unitfapar,
    rd_unitiabs:     rd_unitiabs,
    actnv:           actnv,
    actnv_unitfapar: actnv_unitfapar,
    actnv_unitiabs:  actnv_unitiabs,
    lue:             lue,
  });

  return out_pmodel;
}

function lue_approx( temp, vpd, elv, ca, gstar, ns_star, kmm ){
  ////////////////////////////////////////////////////////////////////
  // Output:   list: 'm' (unitless), 'chi' (unitless)
  // Returns list containing light use efficiency (m) and ci/ci ratio 
  // (chi) based on the approximation of the theoretical relationships
  // of chi with temp, vpd, and elevation. Is now based on SI units as 
  // inputs.
  //------------------------------------------------------------------
  // arguments:
  // temp       deg C, air temperature
  // vpd        Pa, vapour pressure deficit
  // elv        m, elevation above sea level
  // ca         Pa, ambient CO2 partial pressure
  // gstar      Pa, photores. comp. point (Gamma-star)
  // ns_star    (unitless) viscosity correction factor for water
  // kmm        Pa, Michaelis-Menten coeff.

  // Wang-Han Equation:
  var whe = (ee.Image(1.19).add( ee.Image(0.0545).multiply( temp.subtract(25.0) ) ).subtract( ee.Image(0.5).multiply(ee.Image(vpd).log()) ).subtract( ee.Image(8.15e-5).multiply(elv) ) ).exp();

  // leaf-internal-to-ambient CO2 partial pressure (ci/ca) ratio
  var chi = whe.divide( ee.Image(1.0).add(whe) );

  //  m
  var gamma = gstar.divide(ca);
  var m = (chi.subtract(gamma)).divide( chi.add(ee.Image(2).multiply(gamma)) );

  // return output object
  var out_lue = {
    chi : chi,
    m : m,
    n : -9999
  };

  return out_lue;

}


function lue_vpd_c3_simpl( kmm, gstar, ns_star, ca, vpd ){
  ////////////////////////////////////////////////////////////////////
  // Output:   float, ratio of ci/ca (chi)
  // Returns an estimate of leaf internal to ambient CO2
  // partial pressure following the "simple formulation".
  //-----------------------------------------------------------------------
  // arguments
  // kmm     : Pa, Michaelis-Menten coeff.
  // gstar   : Pa, photores. comp. point (Gamma-star)
  // ns_star : (unitless) viscosity correction factor for water
  // ca      : Pa, ambient CO2 partial pressure
  // vpd     : Pa, vapor pressure deficit

  // leaf-internal-to-ambient CO2 partial pressure (ci/ca) ratio
  var xi  = ( (params_gpp.mod(beta)).multiply(kmm).divide(ee.Image(1.6).multiply(ns_star)) ).sqrt();
  var chi = xi.divide(xi.add( vpd.sqrt() ));

  // light use efficiency (m)
  // consistent with this, directly return light-use-efficiency (m)
  var m = ( xi.multiply(ca.add(gstar)).subtract( gstar.multiply(vpd.sqrt()) ) ) .divide ( (xi.multiply(ca.add(ee.Image(2.0).multiply(gstar)))) .add (ee.Image(2.0).multiply(gstar).multiply(vpd.sqrt()) ));

  // n 
  var gamma = gstar.divide(ca);
  var kappa = kmm.divide(ca);
  var n = (chi.add(kappa)).divide(chi.add(ee.Image(2).multiply(gamma)));

  // return output object
  var out_lue = {
    chi : chi,
    m : m,
    n : n
  };
 
  return out_lue;

}

function lue_vpd_c3_full( kmm, gstar, ns_star, ca, vpd ){
  ////////////////////////////////////////////////////////////////////
  // Output:   float, ratio of ci/ca (chi)
  // Features: Returns an estimate of leaf internal to ambient CO2
  //           partial pressure following the "simple formulation".
  //-----------------------------------------------------------------------
  // arguments
  // kmm     : Pa, Michaelis-Menten coeff.
  // gstar   : Pa, photores. comp. point (Gamma-star)
  // ns_star : (unitless) viscosity correction factor for water
  // ca      : Pa, ambient CO2 partial pressure
  // vpd     : Pa, vapor pressure deficit

  // leaf-internal-to-ambient CO2 partial pressure (ci/ca) ratio
  var xi  = ( ( ee.Image(params_pft_gpp.beta).multiply( kmm.add(gstar) ) ).divide( ee.Image(1.6).multiply(ns_star) ) ).sqrt();     // see Eq. 2 in 'Estimation_of_beta.pdf'
  var chi = (gstar.divide(ca)) .add( ( ee.Image(1.0).subtract((gstar).divide(ca) )) .multiply(xi) .divide ( (xi).add(ee.Image(vpd).sqrt()) ));           // see Eq. 1 in 'Estimation_of_beta.pdf'

  
  // Define variable substitutes:
  var vdcg = ca.subtract(gstar);
  var vacg = ca.add(ee.Image(2.0).multiply(gstar));
  var vbkg = ee.Image(params_pft_gpp.beta).multiply(kmm.add(gstar));
  // Check for negatives:
  var vsr = ee.Image(0);
  var m = ee.Image(0);
  var vsr = vsr.where( vbkg.gt(0), ( ee.Image(1.6).multiply(ns_star).multiply(vpd).divide(vbkg) ).sqrt() )
  var m = m.where( vbkg.gt(0), vdcg.divide( vacg.add(ee.Image(3.0).multiply(gstar).multiply(vsr) ) ) )  
  
  // n 
  var gamma = gstar.divide(ca);
  var kappa = kmm.divide(ca);
  var n = (chi.add(kappa)).divide(chi.add(ee.Image(2).multiply(gamma)));
  
  // return output object
  var out_lue = {
    chi : chi,
    m : m,
    n : n
  };

  return out_lue;

}


function lue_c4(){
  ////////////////////////////////////////////////////////////////////
  // Output:   float, ratio of ci/ca (chi)
  // Features: Returns an estimate of leaf internal to ambient CO2
  //           partial pressure following the "simple formulation".
  //-----------------------------------------------------------------------

  // return output object
  var out_lue = {
    chi : ee.Image(-9999),
    m : ee.Image(1),
    n : ee.Image(1)
  };

  return out_lue;

}


function calc_mprime( m ){
  //-----------------------------------------------------------------------
  // Input:  m   (unitless): factor determining LUE
  // Output: mpi (unitless): modiefied m accounting for the co-limitation
  //                         hypothesis after Prentice et al. (2014)
  //-----------------------------------------------------------------------
  var kc = 0.41;          // Jmax cost coefficient

  // square of m-prime (mpi)
  var mprime = m.pow(2).subtract(ee.Image(kc).pow(2.0/3.0).multiply(m.pow(4.0/3.0)) );


  // Check for negatives and take root of square
  var mprime = mprime.where(mprime.gt(0), mprime.sqrt())

  return mprime;
  
}


function co2_to_ca( co2, patm ){
  //-----------------------------------------------------------------------
  // Output:   - ca in units of Pa
  // Features: Converts ca (ambient CO2) from ppm to Pa.
  //-----------------------------------------------------------------------
  // arguments
  // co2  : ambient CO2 in units of ppm
  // patm : monthly atm. pressure, Pa

  var ca = ee.Image(1.e-6).multiply(co2).multiply(patm);         // Pa, atms. CO2
    
  return ca;

}


function ca_to_co2( ca, patm ){
  //-----------------------------------------------------------------------
  // Output:   - co2 in units of Pa
  // Features: Converts ca (ambient CO2) from Pa to ppm.
  //-----------------------------------------------------------------------
  // arguments
  // ca   : ambient CO2 in units of Pa
  // patm : monthly atm. pressure, Pa

  var co2   = ca * ( 1.e6 ) / patm;

  return co2;
  
}


function calc_k( tc, patm ){
  //-----------------------------------------------------------------------
  // Features: Returns the temperature & pressure dependent Michaelis-Menten
  //           coefficient, K (Pa).
  // Ref:      Bernacchi et al. (2001), Improved temperature response 
  //           functions for models of Rubisco-limited photosynthesis, 
  //           Plant, Cell and Environment, 24, 253--259.
  //-----------------------------------------------------------------------
  // arguments
  // tc   : air temperature, deg C 
  // patm : atmospheric pressure, Pa

  // parameters
  var kc25 = 39.97;      // Pa, assuming 25 deg C & 98.716 kPa
  var ko25 = 2.748e4;    // Pa, assuming 25 deg C & 98.716 kPa
  var dhac = 79430;      // J/mol
  var dhao = 36380;      // J/mol
  var kR   = 8.3145;     // J/mol/K
  var kco  = 2.09476e5;  // ppm, US Standard Atmosphere
  var constant2 = ee.Image(298.15).multiply(ee.Image(kR)).multiply((tc.add(ee.Image(273.15)))) ;
  
  var kc = ee.Image(kc25).multiply(( ee.Image(dhac).multiply(tc.subtract(ee.Image(25.0))).divide(constant2) ).exp());
  var ko = ee.Image(ko25).multiply(( ee.Image(dhao).multiply(tc.subtract(ee.Image(25.0))).divide(constant2) ).exp());

  var po = ee.Image(kco).multiply(ee.Image(1e-6)).multiply(patm); // O2 partial pressure
  var k  = kc.multiply( ee.Image(1.0).add(po.divide(ko)) );

  return k;

}


function calc_gstar( tc ){
  //-----------------------------------------------------------------------
  // Features: Returns the temperature-dependent photorespiratory 
  //           compensation point, Gamma star (Pascals), based on constants 
  //           derived from Bernacchi et al. (2001) study. Corresponds
  //           to 'calc_gstar_colin' in pmodel.R.
  // Ref:      Colin's document
  // function return variable
  // real :: gstar   // gamma-star (Pa)
  //-----------------------------------------------------------------------
  // arguments
  // tc  : air temperature (degrees C)

  // local variables
  var gs25 = ee.Image(4.220);    // Pa, assuming 25 deg C & 98.716 kPa)
  var kR   = ee.Image(8.3145);   // J/mol/K
  var dha  = ee.Image(37830);    // J/mol

  // conversion to temperature in Kelvin
  var tk = tc.add(ee.Image(273.15));

  var gstar = gs25.multiply(( dha.divide(kR) ).multiply(( ee.Image(1.0/298.15).subtract(ee.Image(1.0).divide(tk)))).exp());
  
  return gstar;

}


function calc_vcmax25( vcmax, tc ){
  //-----------------------------------------------------------------------
  // Output:   vcmax25  : Vcmax at 25 deg C
  // Features: Returns the temperature-corrected Vcmax at 25 deg C
  // Ref:      Analogue function like 'calc_gstar_gepisat' in Python version
  //           and 'calc_vcmax25_colin' in R version (pmodel.R)
  // function return variable
  // real :: vcmax25  // Vcmax at 25 deg C 
  //-----------------------------------------------------------------------
  // arguments
  // vcmax   // Vcmax at a given temperature tc 
  // tc      // air temperature (degrees C)

  // loal variables
  var dhav = ee.Image(65330);    // J/mol
  var kR   = ee.Image(8.3145);   // J/mol/K

  // conversion to temperature in Kelvin
  var tk = tc.add(ee.Image(273.15));

  var vcmax25 = vcmax.multiply(( ee.Image(-1).multiply(dhav.divide(kR)).multiply( ee.Image(1/298.15).subtract(ee.Image(1).divide(tk)) ) ).exp());

  return vcmax25;
  
}


function calc_patm( elv ){
  //-----------------------------------------------------------------------
  // Features: Returns the atmospheric pressure as a function of elevation
  //           and standard atmosphere (1013.25 hPa)
  // Depends:  - connect_sql
  //           - flux_to_grid
  //           - get_data_point
  //           - get_msvidx
  // Ref:      Allen et al. (1998)
  // function return variable
  // real :: patm    // atmospheric pressure at elevation 'elv', Pa 
  //-----------------------------------------------------------------------
  // argument
  // elv : elevation above sea level, m
    // local variables
  var kPo = 101325;   // standard atmosphere, Pa (Allen, 1973)
  var kTo = 298.15;   // base temperature, K (Prentice, unpublished)
  var kL = 0.0065;    // temperature lapse rate, K/m (Allen, 1973)
  var kG = 9.80665;   // gravitational acceleration, m/s**2 (Allen, 1973)
  var kR = 8.3145;    // universal gas constant, J/mol/K (Allen, 1973)
  var kMa = 0.028963; // molecular weight of dry air, kg/mol (Tsilingiris, 2008)
  var constant1 = (kG*kMa/(kR*kL));
 
  // local variables
  var kPo = ee.Image(101325);   // standard atmosphere, Pa (Allen, 1973)
  var kTo = ee.Image(298.15);   // base temperature, K (Prentice, unpublished)
  var kL = ee.Image(0.0065);    // temperature lapse rate, K/m (Allen, 1973)
  var kG = ee.Image(9.80665);   // gravitational acceleration, m/s**2 (Allen, 1973)
  var kR = ee.Image(8.3145);    // universal gas constant, J/mol/K (Allen, 1973)
  var kMa = ee.Image(0.028963); // molecular weight of dry air, kg/mol (Tsilingiris, 2008)

  // Convert elevation to pressure, Pa:
  var patm = kPo.multiply(ee.Image(1.0).subtract(elv.divide(kTo).multiply(kL)).pow(constant1));
  
  return patm;

}







function calc_viscosity_h2o( tc, patm ){
  //-----------------------------------------------------------------------
  // Features: Calculates viscosity of water at a given temperature and 
  //           pressure.
  // Depends:  density_h2o
  // Ref:      Huber, M. L., R. A. Perkins, A. Laesecke, D. G. Friend, J. V. 
  //           Sengers, M. J. Assael, ..., K. Miyagawa (2009) New 
  //           international formulation for the viscosity of H2O, J. Phys. 
  //           Chem. Ref. Data, Vol. 38(2), pp. 101-125.
  // function return variable
  // real :: viscosity_h2o
  //-----------------------------------------------------------------------
  // arguments
  // tc   : air temperature (tc), degrees C
  // patm : atmospheric pressure (patm), Pa

  // print('----- in calc_viscosity_h2o() ------')
  // print(tc,'tc')
  // print(patm,'patm')

  // local variables
  var tk_ast  = ee.Image(647.096);    // Kelvin
  var rho_ast = ee.Image(322.0);      // kg/m**3
  var mu_ast  = ee.Image(1e-6);       // Pa s

  // Get the density of water, kg/m**3
  var rho = calc_density_h2o( tc, patm );

  // Calculate dimensionless parameters:
  var tbar = ( tc.add(ee.Image(273.15)) ).divide(tk_ast);
  var tbarx = tbar.pow(0.5);
  var tbar2 = tbar.pow(2);
  var tbar3 = tbar.pow(3);
  var rbar = rho.divide(rho_ast);
  // print(rbar,'rbar')

  // Calculate mu0 (Eq. 11 & Table 2, Huber et al., 2009):
  var mu0 = ee.Image(1.67752).add( ee.Image(2.20462).divide(tbar) ).add( ee.Image(0.6366564).divide(tbar2) ).subtract( ee.Image(0.241605).divide(tbar3) );
  var mu0 = ee.Image(1e2).multiply(tbarx).divide(mu0);
  // print(mu0,'mu0')

  // Create Table 3, Huber et al. (2009):
  var h_array = [
                  [ 0.520094, 0.0850895, -1.08374, -0.289555, 0.0, 0.0 ],  // hj0
                  [ 0.222531, 0.999115, 1.88797, 1.26613, 0.0, 0.120573 ], // hj1
                  [ -0.281378, -0.906851, -0.772479, -0.489837, -0.257040, 0.0 ], // hj2
                  [ 0.161913,  0.257399, 0.0, 0.0, 0.0, 0.0 ], // hj3
                  [ -0.0325372, 0.0, 0.0, 0.0698452, 0.0, 0.0 ], // hj4
                  [ 0.0, 0.0, 0.0, 0.0, 0.00872102, 0.0 ], // hj5
                  [ 0.0, 0.0, 0.0, -0.00435673, 0.0, -0.000593264 ] // hj6
                ];

  // Calculate mu1 (Eq. 12 & Table 3, Huber et al., 2009):
  var mu1 = ee.Image(0.0);
  var ctbar = (ee.Image(1.0).divide(tbar)).subtract(ee.Image(1.0));

  for (var i = 0; i < 6; i++) { 
    var coef1 = ctbar.pow(ee.Image(i));
    // print(i,'coef1:',coef1)
    var coef2 = ee.Image(0.0);
    for (var j = 0; j < 7; j++) { 
      // print(j,i,'h_array[j][i]:',h_array[j][i])
      coef2 = coef2.add( ee.Image(h_array[j][i]).multiply((rbar.subtract(ee.Image(1.0))).pow(j)) );
      // print(i,j,'coef2:',coef2)
    }
    mu1 = mu1.add(coef1.multiply(coef2));
    // print(i,'mu1:',mu1)
  }
  mu1 = ( rbar.multiply(mu1) ).exp();
  // print(i, 'mu1', mu1)

  // Calculate mu_bar (Eq. 2, Huber et al., 2009)
  //   assumes mu2 = 1
  var mu_bar = mu0.multiply(mu1);

  // Calculate mu (Eq. 1, Huber et al., 2009)
  var viscosity_h2o = mu_bar.multiply(mu_ast);    // Pa s

  // print('----- END calc_viscosity_h2o() ------')

  return viscosity_h2o;

}



function calc_density_h2o( tc, patm ){
  //-----------------------------------------------------------------------
  // Features: Calculates density of water at a given temperature and 
  //           pressure using the Tumlirz Equation
  // Ref:      F.H. Fisher and O.E Dial, Jr. (1975) Equation of state of 
  //           pure water and sea water, Tech. Rept., Marine Physical 
  //           Laboratory, San Diego, CA.
  // function return variable
  // real :: density_h2o  // density of water, kg/m**3
  //-----------------------------------------------------------------------
  // arguments
  // tc   : air temperature (tc), degrees C
  // patm : atmospheric pressure (patm), Pa

  // Calculate lambda, (bar cm**3)/g:
  var my_lambda = ee.Image(1788.316) .add(
               ee.Image(21.55053).multiply(tc) ).add(
             ee.Image(-0.4695911).multiply(tc).multiply(tc) ).add(
          ee.Image((3.096363e-3)).multiply(tc).multiply(tc).multiply(tc) ).add(
     ee.Image(-1.0*(7.341182e-6)).multiply(tc).multiply(tc).multiply(tc).multiply(tc) );

  // Calculate po, bar
  var po = ee.Image(5918.499) .add( 
              ee.Image(58.05267).multiply(tc) ).add( 
            ee.Image(-1.1253317).multiply(tc).multiply(tc) ).add(
        ee.Image((6.6123869e-3)).multiply(tc).multiply(tc).multiply(tc) ).add(
   ee.Image(-1.0*(1.4661625e-5)).multiply(tc).multiply(tc).multiply(tc).multiply(tc) );

  // Calculate vinf, cm**3/g
  var vinf = ee.Image(0.6980547) .add(
    ee.Image(-1.0*(7.435626e-4)).multiply(tc) ).add(
        ee.Image( (3.704258e-5)).multiply(tc).multiply(tc) ).add(
    ee.Image(-1.0*(6.315724e-7)).multiply(tc).multiply(tc).multiply(tc) ).add(
         ee.Image((9.829576e-9)).multiply(tc).multiply(tc).multiply(tc).multiply(tc) ).add(
   ee.Image(-1.0*(1.197269e-10)).multiply(tc).multiply(tc).multiply(tc).multiply(tc).multiply(tc) ).add(
        ee.Image((1.005461e-12)).multiply(tc).multiply(tc).multiply(tc).multiply(tc).multiply(tc).multiply(tc) ).add(
   ee.Image(-1.0*(5.437898e-15)).multiply(tc).multiply(tc).multiply(tc).multiply(tc).multiply(tc).multiply(tc).multiply(tc) ).add(
        ee.Image( (1.69946e-17)).multiply(tc).multiply(tc).multiply(tc).multiply(tc).multiply(tc).multiply(tc).multiply(tc).multiply(tc) ).add(
   ee.Image(-1.0*(2.295063e-20)).multiply(tc).multiply(tc).multiply(tc).multiply(tc).multiply(tc).multiply(tc).multiply(tc).multiply(tc).multiply(tc) );

  // Convert pressure to bars (1 bar = 100000 Pa)
  var pbar = ee.Image(1e-5).multiply(patm);
  
  // Calculate the specific volume (cm**3 g**-1):
  var vau = vinf.add(my_lambda.divide((po.add(pbar))));

  // Convert to density (g cm**-3) -> 1000 g/kg; 1000000 cm**3/m**3 -> kg/m**3:
  var density_h2o = ee.Image(1e3).divide(vau);

  return density_h2o;

}
