/*

Copyright © 2013, 2017 Battelle Memorial Institute. All Rights Reserved.

NOTICE:  These data were produced by Battelle Memorial Institute (BATTELLE) under Contract No. DE-AC05-76RL01830 with the U.S. Department of Energy (DOE).  For a five year period from May 28, 2013, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, and perform publicly and display publicly, by or on behalf of the Government.  There is provision for the possible extension of the term of this license.  Subsequent to that period or any extension granted, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.  The specific term of the license can be identified by inquiry made to BATTELLE or DOE.  NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR BATTELLE, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY DATA, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*/

/* DO NOT USE FUNCTIONS THAT ARE NOT THREAD SAFE (e.g. rand(), use Util::getModelRand() instead) */

#ifndef __MODEL_DEFINE_H__
#define __MODEL_DEFINE_H__

//#define NUM_JUNCTION_END_TYPES 1 /* we consider only one junction end type */

#include "biocellion.h"

/* define constants to be used inside model functions */
inline REAL radius_from_volume( REAL volume ) {
        return CBRT( volume * 3.0 / ( 4.0 * MY_PI ) );
}
inline REAL volume_agent(REAL radius){
        return (4.0/3.0) * MY_PI * radius*radius*radius;
}
inline REAL surface_agent(REAL radius){
        return 4.0*MY_PI*radius*radius;
}
inline REAL MonodEquation(  REAL Kc , REAL u ) {
        return  u / ( Kc + u ) ;
}


/* MODEL START */
// const char*MODEL_NAME = "ABM-aggregates-box";
const S32 SYSTEM_DIMENSION = 3;
const REAL IF_GRID_SPACING = 500; // also called "h" in the documentation, in micrometers
const REAL EPSILON = 1e-20;

typedef enum _model_rng_type_e {
	MODEL_RNG_UNIFORM,
        MODEL_RNG_UNIFORM_10PERCENT,
        MODEL_RNG_GAUSSIAN,
	NUM_MODEL_RNGS
} model_rng_type_e;

typedef enum _agent_type_e {
	AGENT_MCARRIER,
	AGENT_CELL_A,
	NUM_AGENT_TYPES
} agent_type_e;

typedef enum _agent_mechanic_type_e { 
        MCARRIER_INERT,
        CELL_A_LIVE,
        CELL_A_DEATH,
        NUM_MECH_TYPES,
} agent_mechanic_type_e;

typedef enum _cell_state_real_e {
        CELL_MODEL_REAL_RADIUS,
        CELL_MODEL_REAL_MASS,
        CELL_MODEL_REAL_EPS,
        CELL_MODEL_REAL_UPTAKE_PCT,
        CELL_MODEL_REAL_DX,
        CELL_MODEL_REAL_DY,
        CELL_MODEL_REAL_DZ,
        CELL_MODEL_REAL_STRESS,
        NUM_CELL_STATE_REALS
} cell_state_real_e;

typedef enum _cell_state_int_e { 
        CELL_MODEL_INT_STATE,
        NUM_CELL_STATE_INTS
} cell_state_int_e ;

/* ODE parameters for biomas change */
typedef enum _ode_net_GrowingCell_var_e {
        ODE_NET_VAR_GROWING_CELL_BIOMASS,  // dm/dt = constant
        ODE_NET_VAR_STRESS_TIME,   // dT/dt = heaviside 
        NUM_ODE_NET_VAR_GROWING_CELL
} ode_net_GrowingCell_var_e;

typedef enum _cell_mech_real_e {
        CELL_MECH_REAL_FORCE_X,/* shoving & adhesion */
        CELL_MECH_REAL_FORCE_Y,/* shoving & adhesion */
        CELL_MECH_REAL_FORCE_Z,/* shoving & adhesion */
        CELL_DIVISION_NORMAL_X, // normal to division plane
        CELL_DIVISION_NORMAL_Y, //
        CELL_DIVISION_NORMAL_Z, //
        CELL_MECH_REAL_STRESS,
        NUM_CELL_MECH_REALS
} cell_mech_real_e;

typedef enum _particle_extra_output_real_e {
	PARTICLE_EXTRA_OUTPUT_REAL_RADIUS,
        PARTICLE_EXTRA_OUTPUT_REAL_STRESS,
        PARTICLE_EXTRA_OUTPUT_REAL_ID,
        PARTICEL_EXTRA_OUTPUT_REAL_MicroID,
        PARTICLE_EXTRA_OUTPUT_REAL_VX,
        PARTICLE_EXTRA_OUTPUT_REAL_VY,
        PARTICLE_EXTRA_OUTPUT_REAL_VZ,
	NUM_PARTICLE_EXTRA_OUTPUT_REALS
} particle_extra_output_real_e;

typedef enum _grid_summary_real_e {
        GRID_SUMMARY_REAL_MICROCARRIERS,
        GRID_SUMMARY_REAL_LIVE_CELLS,
        GRID_SUMMARY_REAL_LIVE_ATTACHED,
        GRID_SUMMARY_REAL_DEATH,
        GRID_SUMMARY_REAL_REMOVED,
        GRID_SUMMARY_REAL_MAX_DISP,
        GRID_SUMMARY_REAL_MAX_GROWRATE,
        NUM_GRID_SUMMARY_REALS
} grid_summary_real_e;

typedef enum _junction_end_type_e {
        JUNCTION_END_TYPE_MICROCARRIER,
        JUNCTION_END_TYPE_CELL,
        NUM_JUNCTION_END_TYPES
} junction_end_type_e;

typedef enum _junction_end_real_e {
        JUNCTION_END_REAL_DIR_X,
        JUNCTION_END_REAL_DIR_Y,
        JUNCTION_END_REAL_DIR_Z,
        NUM_JUNCTION_END_REALS
} junction_end_real_e;

typedef enum _grid_model_real_e {
        GRID_MODEL_REAL_DEATHCELL_COUNT,
        NUM_GRID_MODEL_REALS
} grid_model_real_e;


const REAL A_DIFFUSION_COEFF_CELLS[ NUM_AGENT_TYPES ] = { 0.0, 0.0 };
const REAL A_AGENT_FRICIONAL_DRAG[ NUM_AGENT_TYPES ] = { 1.2e-4, 1e-5 }; 
const REAL A_CELL_RADIUS[ NUM_AGENT_TYPES ] = {  30 , 14.5 }; //mean microcarrier radius: 80/2 (from UCLA team); cell radius between 7 and 38 (from UCLA team);
const REAL A_DIVISION_RADIUS[  NUM_AGENT_TYPES ] = { 0.0, 15.0 };
const REAL A_MIN_CELL_RADIUS[ NUM_AGENT_TYPES ] ={ 0.0, 11.0  } ;
const REAL A_MAX_CELL_RADIUS[ NUM_AGENT_TYPES ] ={ 0.0, 16.0  };
const REAL A_MAX_CELL_VOL[NUM_AGENT_TYPES]={ 0.0, 4.0*MY_PI*16.0*16.0*16.0/3.0 };
const REAL A_MIN_CELL_VOL[NUM_AGENT_TYPES]={ 0.0,  523.6 };

// parameters of shoving adhesion
const REAL A_AGENT_SHOVING_SCALE[ NUM_MECH_TYPES ] = {1.0, 1.0, 1.0} ;
const REAL A_AGENT_SHOVING_LIMIT[ NUM_MECH_TYPES ] = { 0.0 , 0.0, 0.0} ;
const REAL A_AGENT_ADHESION_S[NUM_MECH_TYPES][NUM_MECH_TYPES]={{0.0, 0.0, 0.0}, {0.0, 0.0,0.0},{0.0,0.0,0.0}};

// parameters of bonding forces: https://doi.org/10.1371/journal.pone.0191089 
const REAL A_AGENT_STIFFNESS_CM = 0.014739999999999998 ;
const REAL A_AGENT_BOND_S_CM = 0.01 ;
const REAL A_AGENT_BOND_DESTROY_FACTOR_CM = 1.0 ;

const REAL A_AGENT_STIFFNESS_CC = 0.0001 ;
const REAL A_AGENT_BOND_S_CC = 0.20082988502465085 ;
const REAL A_AGENT_BOND_DESTROY_FACTOR_CC = 1.1666666666666667 ;

const REAL A_AGENT_STIFFNESS[NUM_MECH_TYPES][NUM_MECH_TYPES] = {
        {A_AGENT_STIFFNESS_CM,A_AGENT_STIFFNESS_CM,A_AGENT_STIFFNESS_CM},
        {A_AGENT_STIFFNESS_CM,A_AGENT_STIFFNESS_CC,A_AGENT_STIFFNESS_CC},
        {A_AGENT_STIFFNESS_CM,A_AGENT_STIFFNESS_CC,A_AGENT_STIFFNESS_CC}
} ;
const REAL A_AGENT_BOND_S[NUM_MECH_TYPES][NUM_MECH_TYPES] = {
        {A_AGENT_BOND_S_CM,A_AGENT_BOND_S_CM,A_AGENT_BOND_S_CM},
        {A_AGENT_BOND_S_CM,A_AGENT_BOND_S_CC,A_AGENT_BOND_S_CC},
        {A_AGENT_BOND_S_CM,A_AGENT_BOND_S_CC,A_AGENT_BOND_S_CC}
} ;
const REAL A_AGENT_BOND_DESTROY_FACTOR[NUM_MECH_TYPES][NUM_MECH_TYPES] = {
        {A_AGENT_BOND_DESTROY_FACTOR_CM,A_AGENT_BOND_DESTROY_FACTOR_CM,A_AGENT_BOND_DESTROY_FACTOR_CM},
        {A_AGENT_BOND_DESTROY_FACTOR_CM,A_AGENT_BOND_DESTROY_FACTOR_CC,A_AGENT_BOND_DESTROY_FACTOR_CC},
        {A_AGENT_BOND_DESTROY_FACTOR_CM,A_AGENT_BOND_DESTROY_FACTOR_CC,A_AGENT_BOND_DESTROY_FACTOR_CC}
} ;
const REAL A_AGENT_BOND_CREATE_FACTOR[NUM_MECH_TYPES][NUM_MECH_TYPES] = {{0.0,1.1,1.0}, {1.1,1.1,1.0},{1.0,1.0,1.0}};


const REAL A_DENSITY_BIOMASS[ NUM_AGENT_TYPES ] = {8.321827089772318e-17,8.355634512324506e-16 }; // Kg / um^3
const REAL DENSITY_MEDIUM = 1.0e-18   ;// Kg / um^3 
const REAL A_MCARRIER_DENSITY_PER_UB =  1 / (4.0 *4.0 * 4.0 ) ;  // 2000 //0.02 ; // 0.1
const REAL INIT_CELLS_PER_MICROCARRIER = 5; //10;
const REAL A_CELL_D_MAX[ NUM_AGENT_TYPES ] = {150.0*1.25, 20.0*1.25};


const REAL ADHESION_S = 0.01;
const REAL RANDOM_VIBRATION_SCALE = 0.05;

const S32 MECH_INTRCT_ELLIPSOID_MAX_ITERS = 100;
const REAL MECH_INTRCT_ELLIPSOID_EPSILON = 1e-10;

const REAL BASELINE_TIME_STEP_DURATION = 0.00002  ; //0.0001; //  seconds
const REAL STEP_TIME = 1.0 ;
const REAL NUM_STATE_AND_GRID_TIME_STEPS_PER_BASELINE = 1 ;

// cell growth
const REAL DOUBLING_TIME = 0.1 ;   //BASELINE_TIME_STEP_DURATION * 100000 ; // I also used 100000
const REAL ODE_CELL_GROWTH_CONSTANT = A_DENSITY_BIOMASS[1] * ( 4.0 * MY_PI * 15.0*15.0*15.0 / 3.0 ) /( 2.0 * DOUBLING_TIME ); //   4pi/3 ( R_div ^3) * density / (2 * doublingtime)
const REAL STRESS_TRESHOLD = 2.1544346900318822e-07 ;
const REAL STRESS_HILL_EXPONENT = 10.0 ;
const REAL STRESS_TRESHOLD_OVER_N = POW( STRESS_TRESHOLD , STRESS_HILL_EXPONENT ) ;
 
// Bioreactor Geometry
const REAL BIO_RADIUS = 55000 * 0.5; // micrometers
const REAL BIO_HEIGHT = 42895; // micrometers

// Boundary forces 
// U =  eps_B * exp( delta / sigma_B )  for delta >0 ; 0 othersie
// This potential generate de forces  Fx = - dU /dx  
// Assuming microcarriers and cells have the same parameters for simplicity
const REAL EPS_BOUNDARY = 2.0 * 1.0e-9;
const REAL SIG_BOUNDARY = 1.0; 

// drag force 
const REAL DYNAMIC_VISCOSITY =  1.0e-9; //  [micro N][s]/ ( [ micro m] [micro m] ) 

// Verlet integration
const S32 AGENT_TRANSLATION_ROTATION_INTEGRATION_STEPS_PER_BASELINE_TIME_STEP = 1 ;
const REAL AGENT_TRANSLATION_ROTATION_PSEUDO_TIME_STEP_DURATION = BASELINE_TIME_STEP_DURATION / ( (REAL) AGENT_TRANSLATION_ROTATION_INTEGRATION_STEPS_PER_BASELINE_TIME_STEP  ); 

const REAL VELOCITY_DAMPING_TEST  = 1.0; // 1.0; //1.0e-4;

// gravity acceleration
const REAL STANDARD_GRAVITY = 9.80665e-6 ; // micro meter / ( s^2 )

// cell death parameters 
const REAL MECH_STRESS_TRESHOLD_DEATH = 2.1544346900318822e-06 ;
const REAL TIME_TO_DEATH = 0.1 *  DOUBLING_TIME   ;   

/* MODEL END */
#endif/* #ifndef __MODEL_DEFINE_H__ */

