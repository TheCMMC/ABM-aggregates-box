/*

Copyright © 2013, 2017 Battelle Memorial Institute. All Rights Reserved.

NOTICE:  These data were produced by Battelle Memorial Institute (BATTELLE) under Contract No. DE-AC05-76RL01830 with the U.S. Department of Energy (DOE).  For a five year period from May 28, 2013, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, and perform publicly and display publicly, by or on behalf of the Government.  There is provision for the possible extension of the term of this license.  Subsequent to that period or any extension granted, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.  The specific term of the license can be identified by inquiry made to BATTELLE or DOE.  NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR BATTELLE, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY DATA, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*/

/* DO NOT USE FUNCTIONS THAT ARE NOT THREAD SAFE (e.g. rand(), use Util::getModelRand() instead) */

#include "biocellion.h"

#include "model_routine.h"

/* MODEL START */

#include "model_define.h"

extern "C" {
  void cfd_query(double, double, double, double *, double *, double *);
};

/* MODEL END */

using namespace std;

#if HAS_SPAGENT


static void computeAgentTranslation( const VReal& vForce, const VReal& vPos, const JunctionData& junctionData, const MechIntrctData& mechIntrctData, const VReal& FluidVelocity, SpAgentState& state, /* INOUT */ VReal& vDisp ) ;


void ModelRoutine::addSpAgents( const BOOL init, const VIdx& startVIdx, const VIdx& regionVSize, const IfGridBoxData<BOOL>& ifGridHabitableBoxData, Vector<VIdx>& v_spAgentVIdx, Vector<SpAgentState>& v_spAgentState, Vector<VReal>& v_spAgentVOffset ) {/* initialization */
	/* MODEL START */

        REAL xo = REAL ( Info::getDomainSize(0) * IF_GRID_SPACING ) * 0.5 ;
        REAL yo = REAL ( Info::getDomainSize(1) * IF_GRID_SPACING ) * 0.5 ; 

	if( init == true ) {
		const S64 numUBs = regionVSize[0] * regionVSize[1] * regionVSize[2];

		S64 numMcarriers = ( S64 )( ( REAL )numUBs * A_MCARRIER_DENSITY_PER_UB  );
		//S64 numMcarriers = 1 ;

		for( S64 j = 0 ; j < numMcarriers ; j++ ) {

			VReal vPos;
			VIdx vIdx;
			VReal vOffset;
			SpAgentState state;
                        REAL dist = 0.0 ;
                    do {
			for( S32 dim = 0 ; dim < DIMENSION ; dim++ ) {

				REAL randScale = Util::getModelRand( MODEL_RNG_UNIFORM ) ;/* [0.0,1.0) */
				if( randScale >= 1.0 ) {
					randScale = 1.0 - EPSILON;
				}
				CHECK( randScale >= 0.0 );
				CHECK( randScale < 1.0 );
				vPos[dim] = ( REAL )startVIdx[dim] * IF_GRID_SPACING + ( REAL )regionVSize[dim] * IF_GRID_SPACING * randScale;
				//vPos[dim] =  REAL ( Info::getDomainSize( dim ) ) * IF_GRID_SPACING  * 0.5  +   BIO_RADIUS * randScale ;
			}

                        dist = SQRT( (vPos[0] - xo)*(vPos[0] - xo) + (vPos[1] - yo)*(vPos[1] - yo) );
                    } while (  ( dist >= BIO_RADIUS - 290.0 ) || ( vPos[2] >= BIO_HEIGHT - 290.0  ) || ( vPos[2] <= 290.0 ) ) ; 
                             
                        

			Util::changePosFormat1LvTo2Lv( vPos, vIdx, vOffset );


			state.setType( AGENT_MCARRIER );
			state.setModelReal( CELL_MODEL_REAL_RADIUS, A_CELL_RADIUS[AGENT_MCARRIER] );
                        REAL biomass = volume_agent(A_CELL_RADIUS[AGENT_MCARRIER] ) * A_DENSITY_BIOMASS[ AGENT_MCARRIER ] ;
                        state.setModelReal( CELL_MODEL_REAL_MASS, biomass );
                        state.setModelReal( CELL_MODEL_REAL_EPS,  0.0 );
                        state.setModelReal( CELL_MODEL_REAL_UPTAKE_PCT,  1.0 ) ;
                        state.setModelReal( CELL_MODEL_REAL_DX, 0.0 );
                        state.setModelReal( CELL_MODEL_REAL_DY, 0.0 );
                        state.setModelReal( CELL_MODEL_REAL_DZ, 0.0 );
                        state.setModelReal( CELL_MODEL_REAL_STRESS, 0.0 );

                        state.setModelInt( CELL_MODEL_INT_STATE, MCARRIER_INERT ) ; 

			state.setMechIntrctBdrySphere( A_CELL_D_MAX[AGENT_MCARRIER] );

			v_spAgentVIdx.push_back( vIdx );
			v_spAgentState.push_back( state );
			v_spAgentVOffset.push_back( vOffset );


                        // for each microcarrier generate cells
                        S32 numCells =  INIT_CELLS_PER_MICROCARRIER ; // should from Poisson distribution
			for ( S32 i = 0 ; i < numCells ; i++ ) {

                             VReal vPos_c;
                             VIdx vIdx_c;
                             VReal vOffset_c;
                             SpAgentState state_c;

                             REAL randScale = Util::getModelRand( MODEL_RNG_UNIFORM ) ;/* [0.0,1.0) */
                             if( randScale >= 1.0 ) {
                                 randScale = 1.0 - EPSILON;
                             }
                             REAL cellrad = A_MIN_CELL_RADIUS[AGENT_CELL_A] + ( A_CELL_RADIUS[AGENT_CELL_A] - A_MIN_CELL_RADIUS[AGENT_CELL_A] ) * randScale ;
                              
                             REAL rho = A_CELL_RADIUS[AGENT_MCARRIER] + cellrad ;

                             

                             REAL V1, V2, V3, S ;
                             do {
                                 V1 = 2.0 * Util::getModelRand( MODEL_RNG_UNIFORM ) - 1.0 ;
                                 V2 = 2.0 * Util::getModelRand( MODEL_RNG_UNIFORM ) - 1.0 ;
                                 V3 = 2.0 * Util::getModelRand( MODEL_RNG_UNIFORM ) - 1.0 ;

                                 S  = V1*V1 +  V2*V2 + V3*V3 ;
                             }
                             while (  S >= 1.0  || S < 0 ) ;

			     REAL sqrtS = SQRT( S ) ;
			     vPos_c[0] = vPos[0] + rho * V1 / sqrtS  ;
			     vPos_c[1] = vPos[1] + rho * V2 / sqrtS  ;
			     vPos_c[2] = vPos[2] + rho * V3 / sqrtS  ;
 
                             for ( S32 k = 0 ; k < 3; k++ ) { 
                                if  ( vPos_c[k] > 32 * IF_GRID_SPACING )  // 32 ?? change this 
                                    vPos_c[k] =  vPos_c[k] - 32.0 * IF_GRID_SPACING ; 
                                else if (  vPos_c[0] < 0.0 )  
                                    vPos_c[k] =  32.0 * IF_GRID_SPACING - vPos_c[k] ;
 
                             }
			     Util::changePosFormat1LvTo2Lv( vPos_c, vIdx_c, vOffset_c );
                                       

			     state_c.setType( AGENT_CELL_A );
                             state_c.setModelReal( CELL_MODEL_REAL_RADIUS, cellrad );
                             REAL biomass = volume_agent( cellrad )*A_DENSITY_BIOMASS[ AGENT_CELL_A ] ;
                             state_c.setModelReal( CELL_MODEL_REAL_MASS, biomass );
                             state_c.setModelReal( CELL_MODEL_REAL_EPS, 0.0 );
                             state_c.setModelReal( CELL_MODEL_REAL_UPTAKE_PCT, 1.0 ) ;
                             state_c.setModelReal( CELL_MODEL_REAL_DX, 0.0 );
                             state_c.setModelReal( CELL_MODEL_REAL_DY, 0.0 );
                             state_c.setModelReal( CELL_MODEL_REAL_DZ, 0.0 );
                             state_c.setModelReal( CELL_MODEL_REAL_STRESS, 0.0 );

                             state_c.setODEVal(0, ODE_NET_VAR_GROWING_CELL_BIOMASS, biomass );
                             state_c.setODEVal(0, ODE_NET_VAR_STRESS_TIME, 0.0 );

                             state_c.setModelInt( CELL_MODEL_INT_STATE, CELL_A_LIVE ) ; 


                             state_c.setMechIntrctBdrySphere( A_CELL_D_MAX[ AGENT_CELL_A ] );

                             v_spAgentVIdx.push_back( vIdx_c );
                             v_spAgentState.push_back( state_c );
                             v_spAgentVOffset.push_back( vOffset_c );
                             
                            
                        }
		}
	}

	/* MODEL END */

	return;
}

void ModelRoutine::spAgentCRNODERHS( const S32 odeNetIdx, const VIdx& vIdx, const SpAgent& spAgent, const NbrUBEnv& nbrUBEnv, const Vector<double>& v_y, Vector<double>& v_f ) {
	/* MODEL START */

        S32 mtype = spAgent.state.getModelInt( CELL_MODEL_INT_STATE )    ; 
        //cout << " ODE-type " << mtype << endl; 
        if ( mtype == CELL_A_LIVE )  { 
           REAL mech_stress  = spAgent.state.getModelReal( CELL_MODEL_REAL_STRESS  ) ;
           if ( mech_stress < 0.0 )   { 
              mech_stress = -1.0 * mech_stress ;
           }

           //REAL factor = POW( mech_stress , STRESS_HILL_EXPONENT  ); 
           REAL factor = STRESS_TRESHOLD_OVER_N/(STRESS_TRESHOLD_OVER_N + POW( mech_stress , STRESS_HILL_EXPONENT)  );
           REAL r_Growth = ODE_CELL_GROWTH_CONSTANT * factor ; //* biomass;
           v_f[ODE_NET_VAR_GROWING_CELL_BIOMASS]= r_Growth ;

           // cell death: dT/dt = Hevisaide( Stress_i - Stress_0 ) 
           REAL heaviside = 0.0 ;
           if ( mech_stress >  MECH_STRESS_TRESHOLD_DEATH ) {
               heaviside = 1.0;
           }
           v_f[ODE_NET_VAR_STRESS_TIME ]= heaviside ;
        }
        else if ( mtype == CELL_A_DEATH )  { 
           REAL r_Growth = -0.5 *ODE_CELL_GROWTH_CONSTANT  ; //* biomass;
           v_f[ODE_NET_VAR_GROWING_CELL_BIOMASS]= r_Growth ;
           v_f[ODE_NET_VAR_STRESS_TIME ] = 0.0 ;
        }
        else { // error
           v_f[ODE_NET_VAR_GROWING_CELL_BIOMASS]= 0.0 ;
           v_f[ODE_NET_VAR_STRESS_TIME ] = 0.0 ;

        }


	/* MODEL END */

	return;
}

void ModelRoutine::updateSpAgentState( const VIdx& vIdx, const JunctionData& junctionData, const VReal& vOffset, const NbrUBEnv& nbrUBEnv, SpAgentState& state/* INOUT */ ) {
	/* MODEL START */

    //REAL uptakePct = state.getModelReal( CELL_MODEL_REAL_UPTAKE_PCT );
    S32 type = state.getType() ; // id of cell type
     
    //if ( ( type == AGENT_CELL_A ) && (  uptakePct > 0.0 ) ) {
    if ( type == AGENT_CELL_A )  {

        REAL Biomas = state.getODEVal( 0, ODE_NET_VAR_GROWING_CELL_BIOMASS );
        REAL Inert = 0.0;

        REAL cellVol = (Biomas + Inert)/A_DENSITY_BIOMASS[type];

        if ( cellVol > A_MAX_CELL_VOL[type] ) {
            cellVol = A_MAX_CELL_VOL[type] ;
            Biomas = cellVol * A_DENSITY_BIOMASS[type]  - Inert ;
        }
        CHECK( cellVol >= 0.0 );
        REAL newRadius = radius_from_volume( cellVol );

        state.setModelReal(  CELL_MODEL_REAL_RADIUS, newRadius );
        state.setModelReal( CELL_MODEL_REAL_MASS, Biomas );

        //REAL ODE_NET_VAR_STRESS_TIME
        REAL stress_time = state.getODEVal(0, ODE_NET_VAR_STRESS_TIME);
        if ( stress_time >  TIME_TO_DEATH ) {
            // change the cell to cell death type
            state.setModelInt( CELL_MODEL_INT_STATE, CELL_A_DEATH );
            //state.setModelReal( CELL_MODEL_REAL_UPTAKE_PCT, 1.0 ) ;
            state.setODEVal(0, ODE_NET_VAR_STRESS_TIME, 0.0 );
        }

        if ( (newRadius <= A_MIN_CELL_RADIUS[type] ) && (state.getModelReal( CELL_MODEL_REAL_UPTAKE_PCT) == 1.0) ) 
            state.setModelReal( CELL_MODEL_REAL_UPTAKE_PCT, 0.0 ) ;
        else if (( newRadius <= A_MIN_CELL_RADIUS[type] ) && (state.getModelReal( CELL_MODEL_REAL_UPTAKE_PCT) == 0.0))
            state.setModelReal( CELL_MODEL_REAL_UPTAKE_PCT, -1.0 ) ;
    }
 

    /* MODEL END */
     return;
}

void ModelRoutine::spAgentSecretionBySpAgent( const VIdx& vIdx, const JunctionData& junctionData, const VReal& vOffset, const MechIntrctData& mechIntrctData, const NbrUBEnv& nbrUBEnv, SpAgentState& state/* INOUT */, Vector<SpAgentState>& v_spAgentState, Vector<VReal>& v_spAgentVDisp ) {
	/* MODEL START */

	/* nothing to do */

	/* MODEL END */

	return;
}

void ModelRoutine::updateSpAgentBirthDeath( const VIdx& vIdx, const SpAgent& spAgent, const MechIntrctData& mechIntrctData, const NbrUBEnv& nbrUBEnv, BOOL& divide, BOOL& disappear ) {

    /* MODEL START */

    divide = false;
    disappear = false;

    S32 type = spAgent.state.getType() ;

    if ( type == AGENT_CELL_A )  {

        //REAL rnd_num = Util::getModelRand(MODEL_RNG_UNIFORM_10PERCENT); //0.9-1.1
        REAL testrad = A_DIVISION_RADIUS[type] ; //* rnd_num;
        REAL cell_rad = spAgent.state.getModelReal( CELL_MODEL_REAL_RADIUS ) ;

        if ( cell_rad >= testrad ) {
            if ( spAgent.state.getModelInt(CELL_MODEL_INT_STATE) == CELL_A_LIVE ) {     
                for( S32 i = 0 ; i < spAgent.junctionData.getNumJunctions() ; i++ ) {
                     const JunctionEnd& end = spAgent.junctionData.getJunctionEndRef( i );
                     if ( end.getType() == JUNCTION_END_TYPE_MICROCARRIER ) {
                        divide = true;
                        break;
                     }
                }
            }
        }
        else if ( spAgent.state.getModelReal( CELL_MODEL_REAL_UPTAKE_PCT) == -1.0 ) {
            //cout << "seeee " << cell_rad << "  " << type <<  endl;
            disappear = true;
        }
    }

    /* MODEL END */
    return;

}

void ModelRoutine::adjustSpAgent( const VIdx& vIdx, const JunctionData& junctionData, const VReal& vOffset, const MechIntrctData& mechIntrctData, const NbrUBEnv& nbrUBEnv, SpAgentState& state/* INOUT */, VReal& vDisp ) {/* if not dividing or disappearing */

  /* MODEL START */

  VReal vForce ;
  VReal vDragForce ; 
  VReal vFluidV ; 
  VReal vPos =  VReal::ZERO ; 
  
  REAL radius  = state.getModelReal( CELL_MODEL_REAL_RADIUS ); 
  S32 type = state.getType() ;
  REAL xo = REAL ( Info::getDomainSize(0) * IF_GRID_SPACING ) * 0.5 ;
  REAL yo = REAL ( Info::getDomainSize(1) * IF_GRID_SPACING ) * 0.5 ;
  
  REAL Fmag = 0.0 ; 
  
  // retrieve velocities from cfd data based on agent location:
  const int GRID_SIZE = IF_GRID_SPACING;
  double x =  GRID_SIZE * vIdx[0] +  GRID_SIZE*0.5   +  vOffset[0];
  double y =  GRID_SIZE * vIdx[1] +  GRID_SIZE*0.5   +  vOffset[1];
  double z =  GRID_SIZE * vIdx[2] +  GRID_SIZE*0.5   +  vOffset[2];
  double u=0, v=0, w=0;

  cfd_query( (x - xo)*1e-6, (y - yo)*1e-6, z*1e-6, &u, &v, &w); // velocity units m/s
  
  // velocity of fluid in um/s
  vFluidV[0] = u*1e6 * VELOCITY_DAMPING_TEST ;
  vFluidV[1] = v*1e6 * VELOCITY_DAMPING_TEST ;
  vFluidV[2] = w*1e6 * VELOCITY_DAMPING_TEST ;

    
  vForce[0] = mechIntrctData.getModelReal( CELL_MECH_REAL_FORCE_X );
  vForce[1] = mechIntrctData.getModelReal( CELL_MECH_REAL_FORCE_Y );
  vForce[2] = mechIntrctData.getModelReal( CELL_MECH_REAL_FORCE_Z );
  REAL stress = mechIntrctData.getModelReal( CELL_MECH_REAL_STRESS  );    
  //S32 type = state.getType();    

  // Compute force with cylindrical boundary  
  REAL dist = SQRT( (x - xo)*(x-xo) + (y-yo)*(y-yo) );
  REAL delta =  dist + radius - BIO_RADIUS ;
  if ( delta  > 0.0 ) { 
    Fmag =  -( EPS_BOUNDARY / SIG_BOUNDARY ) * EXP( delta / SIG_BOUNDARY ) ; 
    vForce[0] += Fmag * (x - xo) / dist ;
    vForce[1] += Fmag * (y - yo) / dist ;
  }

  // compute forces with bottom and ceiling
  if ( z + radius - BIO_HEIGHT ) { 
    delta = z + radius - BIO_HEIGHT; 
    vForce[2] +=  -( EPS_BOUNDARY / SIG_BOUNDARY ) * EXP( delta / SIG_BOUNDARY );
  }
  else if ( radius - z ) {
    delta = radius ;
    vForce[2] +=  -( EPS_BOUNDARY / SIG_BOUNDARY ) * EXP( delta / SIG_BOUNDARY );
  }
   
  // Gravity and Buoyancy force
  vForce[2] += state.getModelReal(CELL_MODEL_REAL_MASS)*STANDARD_GRAVITY *(DENSITY_MEDIUM/A_DENSITY_BIOMASS[type] - 1.0);    
      
  computeAgentTranslation( vForce, vPos, junctionData, mechIntrctData,  vFluidV,  state /* INOUT */,vDisp ) ;


  // Random movement (Brownian)
  //if ( A_DIFFUSION_COEFF_CELLS[type] > 0.0 ){
  //  REAL F_prw = SQRT( 2*A_DIFFUSION_COEFF_CELLS[type] * dt );
  //  for( S32 dim = 0 ; dim < SYSTEM_DIMENSION ; dim++ )
  //    vDisp[dim]+= F_prw* Util::getModelRand(MODEL_RNG_GAUSSIAN);
  //}

  for( S32 dim = 0 ; dim < DIMENSION ; dim++ ) {/* limit the maximum displacement within a single time step */
    if( FABS( vDisp[dim] ) >= IF_GRID_SPACING ) {
      ERROR( "vDisp[" << dim << "] too large: " << vDisp[dim]  );
    }
  }


  REAL CellVol = volume_agent( radius );
  stress =  ( 0.5 / 3.0 )  * stress / CellVol;

  state.setModelReal( CELL_MODEL_REAL_STRESS, stress );  // update stress
    

  /* MODEL END */

  return;
}

void ModelRoutine::divideSpAgent( const VIdx& vIdx, const JunctionData& junctionData, const VReal& vOffset, const MechIntrctData& mechIntrctData, const NbrUBEnv& nbrUBEnv, SpAgentState& motherState/* INOUT */, VReal& motherVDisp, SpAgentState& daughterState, VReal& daughterVDisp, Vector<BOOL>& v_junctionDivide, BOOL& motherDaughterLinked, JunctionEnd& motherEnd, JunctionEnd& daughterEnd ) {
	/* MODEL START */

    CHECK( ( motherState.getType() == AGENT_CELL_A ) && ( motherState.getModelReal( CELL_MODEL_REAL_UPTAKE_PCT ) > 0.0 )/* live */  );

    REAL radius, radius_dougther  ;
    VReal dir;
    VReal ndir ; // unitary vector that define the division plane
    REAL scale;
    S32 type = motherState.getType();
    //REAL OldVol;
    REAL MotherVol, DougtherVol;
    REAL biomas, mother_biomas, dougther_biomas ;
    REAL inert, mother_inert, dougther_inert;
    REAL rnd_num1 = Util::getModelRand(MODEL_RNG_UNIFORM_10PERCENT);//0.9-1.1

    motherVDisp[0] = mechIntrctData.getModelReal( CELL_MECH_REAL_FORCE_X );
    motherVDisp[1] = mechIntrctData.getModelReal( CELL_MECH_REAL_FORCE_Y );
    motherVDisp[2] = mechIntrctData.getModelReal( CELL_MECH_REAL_FORCE_Z );

    ndir[0] = mechIntrctData.getModelReal( CELL_DIVISION_NORMAL_X );
    ndir[1] = mechIntrctData.getModelReal( CELL_DIVISION_NORMAL_Y );
    ndir[2] = mechIntrctData.getModelReal( CELL_DIVISION_NORMAL_Z );


    biomas = motherState.getModelReal( CELL_MODEL_REAL_MASS );
    inert = motherState.getModelReal( CELL_MODEL_REAL_EPS );
    //OldVol = (biomas + inert)/A_DENSITY_BIOMASS[ type ];

    mother_biomas = 0.5 * biomas *  rnd_num1 ;
    mother_inert = 0.5 * inert * rnd_num1 ;

  
    MotherVol = ( mother_biomas + mother_inert)  / A_DENSITY_BIOMASS[ type ];
    if ( MotherVol < A_MIN_CELL_VOL[type ] ){
       MotherVol = A_MIN_CELL_VOL[type ];
       mother_biomas = MotherVol*A_DENSITY_BIOMASS[type] - mother_inert;
    }

    radius = radius_from_volume( MotherVol );

    //  set the model variable
    motherState.setModelReal( CELL_MODEL_REAL_RADIUS, radius );
    motherState.setModelReal( CELL_MODEL_REAL_MASS, mother_biomas  );
    motherState.setModelReal( CELL_MODEL_REAL_EPS, mother_inert  );
    motherState.setModelReal( CELL_MODEL_REAL_UPTAKE_PCT, 1.0 );
    motherState.setModelInt( CELL_MODEL_INT_STATE, CELL_A_LIVE );

    dougther_biomas = biomas - mother_biomas;
    dougther_inert = inert - mother_inert;

    DougtherVol = (dougther_biomas + dougther_inert)/A_DENSITY_BIOMASS[type ];
    if ( DougtherVol < A_MIN_CELL_VOL[type] ){
       DougtherVol = A_MIN_CELL_VOL[type] ;
       dougther_biomas = DougtherVol*A_DENSITY_BIOMASS[type];
    }

    radius_dougther = radius_from_volume( DougtherVol );
    // set the model variable
    daughterState.setType( type );
    daughterState.setModelReal( CELL_MODEL_REAL_RADIUS, radius_dougther ) ; 
    daughterState.setModelReal( CELL_MODEL_REAL_MASS, dougther_biomas );
    daughterState.setModelReal( CELL_MODEL_REAL_EPS, dougther_inert );
    daughterState.setModelReal( CELL_MODEL_REAL_UPTAKE_PCT, 1.0 );
    daughterState.setModelInt( CELL_MODEL_INT_STATE, CELL_A_LIVE );

    // change the new biomass on the  ODEs
    motherState.setODEVal(0, ODE_NET_VAR_GROWING_CELL_BIOMASS, mother_biomas);
    daughterState.setODEVal(0, ODE_NET_VAR_GROWING_CELL_BIOMASS, dougther_biomas);
    daughterState.setODEVal(0, ODE_NET_VAR_STRESS_TIME, 0.0 );

    // mechanical interaction range
    motherState.setMechIntrctBdrySphere( A_CELL_D_MAX[ AGENT_CELL_A ] );
    daughterState.setMechIntrctBdrySphere( A_CELL_D_MAX[ AGENT_CELL_A ] );
    

    // divide in a random direction, other ways are possibl
    VReal rdir = VReal::ZERO;  // random vector in 3D

    do {
        scale = 0.0;
        for( S32 dim = 0; dim < SYSTEM_DIMENSION; dim++ ) {
            rdir[dim] = Util::getModelRand( MODEL_RNG_UNIFORM ) - 1.0;
            scale += rdir[dim] * rdir[dim];
        }
    } while( scale > 1.0 );

    REAL rdir_dot_ndir = ndir[0]*rdir[0] + ndir[1]*rdir[1] + ndir[2]*rdir[2] ;
       
    dir =  rdir - ( ndir * rdir_dot_ndir ) ;  
    scale = dir.length() ;
    dir = dir / scale ; 


    //radius = 0.5* A_AGENT_SHOVING_SCALE[ type_id] *(daughterState.getRadius()+motherState.getRadius());
    radius = 0.25*(daughterState.getModelReal(CELL_MODEL_REAL_RADIUS) + motherState.getModelReal(CELL_MODEL_REAL_RADIUS));
    motherVDisp += dir * radius;
    daughterVDisp -= dir * radius;


    for( S32 dim = 0 ; dim < SYSTEM_DIMENSION ; dim++ ) {/* limit the maximum displacement  */
        if( motherVDisp[dim] > A_MAX_CELL_RADIUS[type] )
            motherVDisp[dim] = A_MAX_CELL_RADIUS[type] ;
        else if( motherVDisp[dim] < ( A_MAX_CELL_RADIUS[type] * -1.0 ) )
            motherVDisp[dim] = A_MAX_CELL_RADIUS[type] * -1.0;

        if( daughterVDisp[dim] > A_MAX_CELL_RADIUS[type]  )
            daughterVDisp[dim] = A_MAX_CELL_RADIUS[type] ;
        else if( daughterVDisp[dim] < ( A_MAX_CELL_RADIUS[type]  * -1.0 ) )
            daughterVDisp[dim] = A_MAX_CELL_RADIUS[type] * -1.0;
    }

    motherDaughterLinked = true;
    motherEnd.setType( JUNCTION_END_TYPE_CELL ); 
    daughterEnd.setType( JUNCTION_END_TYPE_CELL ); 

    /* MODEL END */

    return;
}


// Verlet/leapfrog algorithm for agents translation
static void computeAgentTranslation( const VReal& vForce, const VReal& vPos, const JunctionData& junctionData, const MechIntrctData& mechIntrctData, const VReal& FluidVelocity,  SpAgentState& state, /* INOUT */ VReal& vDisp ) {

  REAL m = state.getModelReal( CELL_MODEL_REAL_MASS );
  S32 type = state.getType() ;
  REAL radius  = state.getModelReal( CELL_MODEL_REAL_RADIUS ); 
    
  VReal oldStaggeredVLinear;
  VReal newStaggeredVLinear;
  VReal vLinear;
  VReal oldVDisp;
  VReal newVDisp;
  VReal tmpVForce = vForce ;


  if ( Info::getCurBaselineTimeStep() > 0 ) {
    oldStaggeredVLinear[0] = state.getModelReal( CELL_MODEL_REAL_DX ) ; // um/s
    oldStaggeredVLinear[1] = state.getModelReal( CELL_MODEL_REAL_DY ) ; // um/s
    oldStaggeredVLinear[2] = state.getModelReal( CELL_MODEL_REAL_DZ ) ; // um/s
  }
  else {
    
    oldStaggeredVLinear = FluidVelocity  -  ( vForce ) * AGENT_TRANSLATION_ROTATION_PSEUDO_TIME_STEP_DURATION / ( m * 2.0 ) ;
    //oldStaggeredVLinear = VReal::ZERO;
    //cout << " m " << m;
    //cout << " Fluid " << FluidVelocity;
    //cout << " vForce " << vForce ;
    //cout << " timestep  " << AGENT_TRANSLATION_ROTATION_PSEUDO_TIME_STEP_DURATION  ;
    //cout << " innertia " <<  AGENT_TRANSLATION_ROTATION_PSEUDO_TIME_STEP_DURATION / ( m * 2.0 )  ;
    //cout << endl ;
     
  }

  vDisp = VReal::ZERO;
  oldVDisp = oldStaggeredVLinear * AGENT_TRANSLATION_ROTATION_PSEUDO_TIME_STEP_DURATION;/* v_{n+1/2} * deltaT */


  for( S32 i = 0 ; i < AGENT_TRANSLATION_ROTATION_INTEGRATION_STEPS_PER_BASELINE_TIME_STEP ; i++ ) { 

    VReal G = ( FluidVelocity - oldStaggeredVLinear ) * DYNAMIC_VISCOSITY * 6.0 * MY_PI * radius  ; // drag force
    newStaggeredVLinear = oldStaggeredVLinear + ( tmpVForce + G ) * ( AGENT_TRANSLATION_ROTATION_PSEUDO_TIME_STEP_DURATION / m );/* v_{n+1/2} = v_{n-1/2} + ( F_n - G_{n-1/2}) * (deltaT/m), G: damping */
  
    newVDisp = newStaggeredVLinear * AGENT_TRANSLATION_ROTATION_PSEUDO_TIME_STEP_DURATION;/* v_{n+1/2} * deltaT */

    vLinear = ( newVDisp + oldVDisp ) / ( AGENT_TRANSLATION_ROTATION_PSEUDO_TIME_STEP_DURATION * 2.0 );/* v_n = (pos_{n+1} - pos_{n-1}) / (deltaT * 2.0) */

    G =    ( FluidVelocity - vLinear ) * DYNAMIC_VISCOSITY * 6.0 * MY_PI * radius  ; 
    newStaggeredVLinear = oldStaggeredVLinear + ( tmpVForce + G ) * ( AGENT_TRANSLATION_ROTATION_PSEUDO_TIME_STEP_DURATION / m );/* v_{n+1/2} = v_{n-1/2} + ( F_n - G_n) * (deltaT/m) */

    newVDisp = newStaggeredVLinear * AGENT_TRANSLATION_ROTATION_PSEUDO_TIME_STEP_DURATION;/* v_{n+1/2} * deltaT */

    if ( type == -1 ) {
    cout << " -------------------" << endl;
    cout << " type  " <<  type << endl ;
    cout << " FluidVelocity " <<  FluidVelocity  << endl;
    cout << " delta V  " << FluidVelocity - oldStaggeredVLinear << endl;
    cout << " Force " << tmpVForce  << endl;
    cout << " G_Half " << ( FluidVelocity - oldStaggeredVLinear ) * DYNAMIC_VISCOSITY * 6.0 * MY_PI * radius << endl ;
    cout << " G " <<  G   << endl;
    cout << " acceler " << ( tmpVForce + G ) * ( AGENT_TRANSLATION_ROTATION_PSEUDO_TIME_STEP_DURATION / m ) << endl;
    cout << " Old Velocity " << oldStaggeredVLinear << endl;
    cout << " vLinear " << vLinear << endl ;
    cout << " New Velocity " <<  newStaggeredVLinear << endl;
    cout << " newDisp " << newVDisp << endl;
    }

    oldStaggeredVLinear = newStaggeredVLinear;
    oldVDisp = newVDisp;
    vDisp += newVDisp;

  }
  
 
  state.setModelReal( CELL_MODEL_REAL_DX, newStaggeredVLinear[0] );  // displacement
  state.setModelReal( CELL_MODEL_REAL_DY, newStaggeredVLinear[1] );
  state.setModelReal( CELL_MODEL_REAL_DZ, newStaggeredVLinear[2] );

}
#endif

