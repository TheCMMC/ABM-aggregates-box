/*

Copyright © 2013, 2017 Battelle Memorial Institute. All Rights Reserved.

NOTICE:  These data were produced by Battelle Memorial Institute (BATTELLE) under Contract No. DE-AC05-76RL01830 with the U.S. Department of Energy (DOE).  For a five year period from May 28, 2013, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, and perform publicly and display publicly, by or on behalf of the Government.  There is provision for the possible extension of the term of this license.  Subsequent to that period or any extension granted, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.  The specific term of the license can be identified by inquiry made to BATTELLE or DOE.  NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR BATTELLE, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY DATA, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*/

/* DO NOT USE FUNCTIONS THAT ARE NOT THREAD SAFE (e.g. rand(), use Util::getModelRand() instead) */

#include "biocellion.h"

#include "model_routine.h"

/* MODEL START */

#include "model_define.h"

// extern "C" {
//   void cfd_query(double, double, double, double *, double *, double *);
// };

/* MODEL END */

using namespace std;

#if HAS_SPAGENT


void ModelRoutine::addSpAgents( const BOOL init, const VIdx& startVIdx, const VIdx& regionVSize, const IfGridBoxData<BOOL>& ifGridHabitableBoxData, Vector<VIdx>& v_spAgentVIdx, Vector<SpAgentState>& v_spAgentState, Vector<VReal>& v_spAgentVOffset ) {/* initialization */
    /* MODEL START */

    if( init == true ) {
        const S64 numUBs = regionVSize[0] * regionVSize[1] * regionVSize[2];
        // S64 numMcarriers = ( S64 )( ( REAL )numUBs * A_MCARRIER_DENSITY_PER_UB  );
        S64 numMcarriers = 1 ;

        for( S64 j = 0 ; j < numMcarriers ; j++ ) {
            VReal vPos;
            VIdx vIdx;
            VReal vOffset;
            SpAgentState state;
        
            for( S32 dim = 0 ; dim < DIMENSION ; dim++ ) {
                REAL randScale = Util::getModelRand( MODEL_RNG_UNIFORM ) ;/* [0.0,1.0) */
                if( randScale >= 1.0 ) {
                    randScale = 1.0 - EPSILON;
                }
                CHECK( randScale >= 0.0 );
                CHECK( randScale < 1.0 );
                vPos[dim] = ( REAL )startVIdx[dim] * IF_GRID_SPACING + ( REAL )regionVSize[dim] * IF_GRID_SPACING * randScale;
            }      
        
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
            S32 numCells =  INIT_CELLS_PER_MICROCARRIER ; // should form Poisson distribution
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

	/* MODEL END */

	return;
}

void ModelRoutine::updateSpAgentState( const VIdx& vIdx, const JunctionData& junctionData, const VReal& vOffset, const NbrUBEnv& nbrUBEnv, SpAgentState& state/* INOUT */ ) {
	/* MODEL START */

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

    /* MODEL END */
    return;

}

void ModelRoutine::adjustSpAgent( const VIdx& vIdx, const JunctionData& junctionData, const VReal& vOffset, const MechIntrctData& mechIntrctData, const NbrUBEnv& nbrUBEnv, SpAgentState& state/* INOUT */, VReal& vDisp ) {/* if not dividing or disappearing */

  /* MODEL START */


  // Random movement (Brownian)
  //if ( A_DIFFUSION_COEFF_CELLS[type] > 0.0 ){
  //  REAL F_prw = SQRT( 2*A_DIFFUSION_COEFF_CELLS[type] * dt );
  //  for( S32 dim = 0 ; dim < SYSTEM_DIMENSION ; dim++ )
  //    vDisp[dim]+= F_prw* Util::getModelRand(MODEL_RNG_GAUSSIAN);
  //}


  /* MODEL END */

  return;
}

void ModelRoutine::divideSpAgent( const VIdx& vIdx, const JunctionData& junctionData, const VReal& vOffset, const MechIntrctData& mechIntrctData, const NbrUBEnv& nbrUBEnv, SpAgentState& motherState/* INOUT */, VReal& motherVDisp, SpAgentState& daughterState, VReal& daughterVDisp, Vector<BOOL>& v_junctionDivide, BOOL& motherDaughterLinked, JunctionEnd& motherEnd, JunctionEnd& daughterEnd ) {
	/* MODEL START */

    /* MODEL END */

    return;
}

#endif

