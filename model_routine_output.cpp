/*

Copyright © 2013, 2017 Battelle Memorial Institute. All Rights Reserved.

NOTICE:  These data were produced by Battelle Memorial Institute (BATTELLE) under Contract No. DE-AC05-76RL01830 with the U.S. Department of Energy (DOE).  For a five year period from May 28, 2013, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, and perform publicly and display publicly, by or on behalf of the Government.  There is provision for the possible extension of the term of this license.  Subsequent to that period or any extension granted, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.  The specific term of the license can be identified by inquiry made to BATTELLE or DOE.  NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR BATTELLE, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY DATA, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*/

/* DO NOT USE FUNCTIONS THAT ARE NOT THREAD SAFE (e.g. rand(), use Util::getModelRand() instead) */

#include "biocellion.h"

#include "model_routine.h"

/* UESR START */

#include "model_define.h"

/* UESR END */

using namespace std;

#if HAS_SPAGENT
void ModelRoutine::updateSpAgentOutput( const VIdx& vIdx, const SpAgent& spAgent, REAL& color, Vector<REAL>& v_extraReal, Vector<VReal>& v_extraVReal ) {
	/* MODEL START */

	color = spAgent.state.getModelInt( CELL_MODEL_INT_STATE );
	CHECK( v_extraReal.size() == NUM_PARTICLE_EXTRA_OUTPUT_REALS );
        REAL radius = spAgent.state.getModelReal( CELL_MODEL_REAL_RADIUS );   
	v_extraReal[PARTICLE_EXTRA_OUTPUT_REAL_RADIUS] = radius;

        // print the Pressure
        //REAL CellVol = volume_agent( radius );
        REAL stress = spAgent.state.getModelReal( CELL_MODEL_REAL_STRESS );

        // print the id of supporting microcarrier, -1 for microcarroers and -1 detached cells
        S64 microID = -1 ;
        if ( spAgent.state.getType() == AGENT_CELL_A ) { 
           for ( S32 i = 0 ; i <  spAgent.junctionData.getNumJunctions(); i++ ) {
              JunctionEnd end = spAgent.junctionData.getJunctionEndRef(i);
              if ( end.getType() == JUNCTION_END_TYPE_MICROCARRIER  ) {
                 microID = spAgent.junctionData.getOtherEndId( i ) ;
              }
           }
        }
        v_extraReal[ PARTICEL_EXTRA_OUTPUT_REAL_MicroID ] = REAL( microID ) ;
        v_extraReal[ PARTICLE_EXTRA_OUTPUT_REAL_STRESS ] =  stress;
        v_extraReal[PARTICLE_EXTRA_OUTPUT_REAL_ID] = REAL( spAgent.junctionData.getCurId() )  ;
        v_extraReal[PARTICLE_EXTRA_OUTPUT_REAL_VX] = spAgent.state.getModelReal( CELL_MODEL_REAL_DX ) / BASELINE_TIME_STEP_DURATION ;
        v_extraReal[PARTICLE_EXTRA_OUTPUT_REAL_VY] = spAgent.state.getModelReal( CELL_MODEL_REAL_DY ) / BASELINE_TIME_STEP_DURATION ;
        v_extraReal[PARTICLE_EXTRA_OUTPUT_REAL_VZ] = spAgent.state.getModelReal( CELL_MODEL_REAL_DZ ) / BASELINE_TIME_STEP_DURATION ;


	/* MODEL END */

	return;
}
#endif

void ModelRoutine::updateSummaryVar( const VIdx& vIdx, const NbrUBAgentData& nbrUBAgentData, const NbrUBEnv& nbrUBEnv, Vector<REAL>& v_realVal/* [elemIdx] */, Vector<S32>& v_intVal/* [elemIdx] */ ) {
   /* MODEL START */

   CHECK( v_realVal.size() == NUM_GRID_SUMMARY_REALS );
   CHECK( v_intVal.size() == 0 );

   const UBAgentData& ubAgentData = *( nbrUBAgentData.getConstPtr( 0, 0, 0 ) );

   REAL count_micro = 0.0;
   REAL count_live = 0.0;
   REAL count_attached = 0.0;
   REAL count_death = 0.0;
   REAL max_disp = 0.0;
   REAL max_fact = -1.0 ;

   /* Count the number of cells placed in the Simulation Domain */
   for (S32 i = 0 ; i < ( S32 )ubAgentData.v_spAgent.size() ; i++ ) {
          
      S32 mtype = ubAgentData.v_spAgent[i].state.getModelInt( CELL_MODEL_INT_STATE );

      if ( mtype == MCARRIER_INERT ) 
         count_micro += 1.0 ;
         if ( ubAgentData.v_spAgent[i].state.getType() == AGENT_MCARRIER ) {
            count_attached  +=  ( REAL) ubAgentData.v_spAgent[i].junctionData.getNumJunctions() ; 
         }
      else if ( mtype == CELL_A_LIVE ) { 
         count_live += 1.0;
      }
      else if ( mtype == CELL_A_DEATH )
         count_death += 1.0;

      REAL dx = FABS(  ubAgentData.v_spAgent[i].state.getModelReal( CELL_MODEL_REAL_DX )) ;
      if ( dx > max_disp )  max_disp = dx ;
          
      REAL dy = FABS(  ubAgentData.v_spAgent[i].state.getModelReal( CELL_MODEL_REAL_DY )) ;
      if ( dy > max_disp )  max_disp = dy ;

      REAL dz = FABS(  ubAgentData.v_spAgent[i].state.getModelReal( CELL_MODEL_REAL_DZ )) ;
      if ( dz > max_disp )  max_disp = dz ;

      REAL mech_stress  = ubAgentData.v_spAgent[i].state.getModelReal( CELL_MODEL_REAL_STRESS  ) ;
      REAL factor =  STRESS_TRESHOLD*STRESS_TRESHOLD / ( STRESS_TRESHOLD*STRESS_TRESHOLD +  mech_stress*mech_stress  );
      if ( factor > max_fact ) max_fact = factor ;

   }


 
   /* URID_SUMMARY_REAL_LIVE_CELLS is set in model_routine_config.cpp */
   v_realVal[ GRID_SUMMARY_REAL_MICROCARRIERS ] = count_micro;
   v_realVal[ GRID_SUMMARY_REAL_LIVE_CELLS ] = count_live;
   v_realVal[ GRID_SUMMARY_REAL_LIVE_ATTACHED ] = count_attached ;
   v_realVal[ GRID_SUMMARY_REAL_DEATH ] = count_death;
   v_realVal[ GRID_SUMMARY_REAL_REMOVED ] = nbrUBEnv.getModelReal(0,0,0,GRID_MODEL_REAL_DEATHCELL_COUNT )  ;
   v_realVal[GRID_SUMMARY_REAL_MAX_DISP ] = max_disp * BASELINE_TIME_STEP_DURATION ;
   v_realVal[GRID_SUMMARY_REAL_MAX_GROWRATE ]  = max_fact ; 
   /* MODEL END */


   /* MODEL END */
   return;
}

