/*

Copyright © 2013, 2017 Battelle Memorial Institute. All Rights Reserved.

NOTICE:  These data were produced by Battelle Memorial Institute (BATTELLE) under Contract No. DE-AC05-76RL01830 with the U.S. Department of Energy (DOE).  For a five year period from May 28, 2013, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, and perform publicly and display publicly, by or on behalf of the Government.  There is provision for the possible extension of the term of this license.  Subsequent to that period or any extension granted, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.  The specific term of the license can be identified by inquiry made to BATTELLE or DOE.  NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR BATTELLE, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY DATA, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*/

/* DO NOT USE FUNCTIONS THAT ARE NOT THREAD SAFE (e.g. rand(), use Util::getModelRand() instead) */

#include "biocellion.h"

#include "model_routine.h"

#include "model_define.h"

using namespace std;

void ModelRoutine::initIfGridVar( const VIdx& vIdx, const UBAgentData& ubAgentData, UBEnv& ubEnv ) {
	/* MODEL START */
	for( S32 sol = 0 ; sol < NUM_GRID_MODEL_REALS ; sol++ ) {
		ubEnv.setModelReal( sol, 0.0 );
	}
	/* MODEL END */

	return;
}

void ModelRoutine::initIfSubgridKappa( const S32 pdeIdx, const VIdx& vIdx, const VIdx& ifSubgridVOffset, const UBAgentData& ubAgentData, const UBEnv& ubEnv, REAL& kappa ) {/* relevant only if v_phiOutputDivideByKappa[pdeIdx] is set to true in updateFileOutputInfo() */
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridVar( const BOOL pre, const S32 iter, const VIdx& vIdx, const NbrUBAgentData& nbrUBAgentData, NbrUBEnv& nbrUBEnv/* [INOUT] */ ) {
	/* MODEL START */

   if ( pre == true ) {
      REAL gridcount =  nbrUBEnv.getModelReal(0,0,0,GRID_MODEL_REAL_DEATHCELL_COUNT ) ;  
      const UBAgentData& ubAgentData = *( nbrUBAgentData.getConstPtr(0,0,0) ) ; 

      for( S32 i = 0 ; i < ( S32 )ubAgentData.v_spAgent.size() ; i++ ) {
         const SpAgent& spAgent = ubAgentData.v_spAgent[i];
         agentType_t type = spAgent.state.getType();

         if ( type == AGENT_CELL_A )  {
            REAL cell_rad = spAgent.state.getModelReal( CELL_MODEL_REAL_RADIUS ) ;

            if ( cell_rad <= A_MIN_CELL_RADIUS[type] ) {
               gridcount += 1.0 ;
            }
         }
      }

      nbrUBEnv.setModelReal(0,0,0,GRID_MODEL_REAL_DEATHCELL_COUNT, gridcount) ;
   }
   /* MODEL END */

   return;
}

void ModelRoutine::updateIfSubgridKappa( const S32 pdeIdx, const VIdx& vIdx, const VIdx& ifSubgridVOffset, const UBAgentData& ubAgentData, const UBEnv& ubEnv, REAL& kappa ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfSubgridAlpha( const S32 elemIdx, const VIdx& vIdx, const VIdx& ifSubgridVOffset, const UBAgentData& ubAgentData, const UBEnv& ubEnv, REAL& alpha/* decay (-) */ ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfSubgridBetaInIfRegion( const S32 elemIdx, const S32 dim, const VIdx& vIdx0, const VIdx& ifSubgridVOffset0, const UBAgentData& ubAgentData0, const UBEnv& ubEnv0, const VIdx& vIdx1, const VIdx& ifSubgridVOffset1, const UBAgentData& ubAgentData1, const UBEnv& ubEnv1, REAL& beta ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfSubgridBetaPDEBufferBdry( const S32 elemIdx, const S32 dim, const VIdx& vIdx, const VIdx& ifSubgridVOffset, const UBAgentData& ubAgentData, const UBEnv& ubEnv, REAL& beta ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfSubgridBetaDomainBdry( const S32 elemIdx, const S32 dim, const VIdx& vIdx, const VIdx& ifSubgridVOffset, const UBAgentData& ubAgentData, const UBEnv& ubEnv, REAL& beta ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

#if 1/* HAS_FLOW */
void ModelRoutine::updateIfSubgridAdvVelInIfRegion( const S32 pdeIdx, const S32 dim, const VIdx& vIdx0, const VIdx& ifSubgridVOffset0, const UBAgentData& ubAgentData0, const UBEnv& ubEnv0, const VIdx& vIdx1, const VIdx& ifSubgridVOffset1, const UBAgentData& ubAgentData1, const UBEnv& ubEnv1, REAL& advVel ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfSubgridAdvVelPDEBufferBdry( const S32 pdeIdx, const S32 dim, const VIdx& vIdx, const VIdx& ifSubgridVOffset, const UBAgentData& ubAgentData, const UBEnv& ubEnv, REAL& advVel ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfSubgridAdvVelDomainBdry( const S32 pdeIdx, const S32 dim, const VIdx& vIdx, const VIdx& ifSubgridVOffset, const UBAgentData& ubAgentData, const UBEnv& ubEnv, REAL& advVel ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}
#endif

void ModelRoutine::updateIfSubgridRHSLinear( const S32 elemIdx, const VIdx& vIdx, const VIdx& ifSubgridVOffset, const UBAgentData& ubAgentData, const UBEnv& ubEnv, REAL& rhs/* uptake(-) and secretion (+) */ ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfSubgridKappaBCVal( const S32 pdeIdx, const S32 dim, const VIdx& vIdx, const VIdx& ifSubgridVOffset, const UBAgentData& ubAgentData, const UBEnv& ubEnv, REAL& kappa ) {/* for boundary condition */
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfSubgridDirichletBCVal( const S32 elemIdx, const S32 dim, const VIdx& vIdx, const VIdx& ifSubgridVOffset, const UBAgentData& ubAgentData, const UBEnv& ubEnv, REAL& bdryPhi ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfSubgridNeumannBCVal( const S32 elemIdx, const S32 dim, const VIdx& vIdx, const VIdx& ifSubgridVOffset, const UBAgentData& ubAgentData, const UBEnv& ubEnv, REAL& bdrySlope ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfSubgridDirichletOutflowBCVal( const S32 elemIdx, const S32 dim, const VIdx& vIdx, const VIdx& ifSubgridVOffset, const UBAgentData& ubAgentData, const UBEnv& ubEnv, REAL& bdryPhi ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfSubgridNeumannOutflowBCVal( const S32 elemIdx, const S32 dim, const VIdx& vIdx, const VIdx& ifSubgridVOffset, const UBAgentData& ubAgentData, const UBEnv& ubEnv, REAL& bdrySlope ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::adjustIfSubgridRHSTimeDependentLinear( const S32 elemIdx, const VIdx& vIdx, const VIdx& ifSubgridVOffset, const UBAgentData& ubAgentData, const UBEnvModelVar& ubEnvModelVar, const REAL phi, REAL& rhs/* INOUT */ ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfSubgridRHSTimeDependentSplitting( const S32 pdeIdx, const VIdx& vIdx, const VIdx& ifSubgridVOffset, const UBAgentData& ubAgentData, const UBEnvModelVar& ubEnvModelVar, const Vector<double>& v_phi/* [idx] */, Vector<double>& v_rhs/* [idx], uptake(-) and secretion (+) */ ) {/* for Wnt & SFRP */
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updateIfGridAMRTags( const VIdx& vIdx, const NbrUBAgentData& nbrUBAgentData, const NbrUBEnv& nbrUBEnv, Vector<S32>& v_finestLevel/* [pdeIdx] */ ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::initPDEBufferGridPhi( const S32 pdeIdx, const VIdx& startVIdx, const VIdx& pdeBufferBoxVSize, Vector<REAL>& v_phi/* [idx] */ ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::initPDEBufferGridKappa( const S32 pdeIdx, const VIdx& startVIdx, const VIdx& pdeBufferBoxVSize, REAL& kappa ) {/* relevant only if v_phiOutputDivideByKappa[pdeIdx] is set to true in updateFileOutputInfo() */
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferGridKappa( const S32 pdeIdx, const VIdx& startVIdx, const VIdx& pdeBufferBoxVSize, REAL& kappa ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferGridAlpha( const S32 elemIdx, const VIdx& startVIdx, const VIdx& pdeBufferBoxVSize, REAL& alpha/* decay (-) */ ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferGridBetaInPDEBufferRegion( const S32 elemIdx, const S32 dim, const VIdx& startVIdx0, const VIdx& startVIdx1, const VIdx& pdeBufferBoxVSize, REAL& beta ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferGridBetaDomainBdry( const S32 elemIdx, const S32 dim, const VIdx& startVIdx, const VIdx& pdeBufferBoxVSize, REAL& beta ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

#if 1/* HAS_FLOW */
void ModelRoutine::updatePDEBufferGridAdvVelInPDEBufferRegion( const S32 pdeIdx, const S32 dim, const VIdx& startVIdx0, const VIdx& startVIdx1, const VIdx& pdeBufferBoxVSize, REAL& advVel ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferGridAdvVelDomainBdry( const S32 pdeIdx, const S32 dim, const VIdx& startVIdx, const VIdx& pdeBufferBoxVSize, REAL& advVel ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}
#endif

void ModelRoutine::updatePDEBufferGridRHSLinear( const S32 elemIdx, const VIdx& startVIdx, const VIdx& pdeBufferBoxVSize, const REAL phi, REAL& rhs/* uptake(-) and secretion (+) */ ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferGridKappaBCVal( const S32 pdeIdx, const S32 dim, const VIdx& startVIdx, const VIdx& pdeBufferBoxVSize, REAL& kappa ) {/* for boundary condition */
	/* MODEL START */

	ERROR( "unimplmented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferGridDirichletBCVal( const S32 elemIdx, const S32 dim, const VIdx& startVIdx, const VIdx& pdeBufferFaceVSize, REAL& bdryPhi ) {
	/* MODEL START */

	ERROR( "unimplmented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferGridNeumannBCVal( const S32 elemIdx, const S32 dim, const VIdx& startVIdx, const VIdx& pdeBufferFaceVSize, REAL& bdrySlope ) {
	/* MODEL START */

	ERROR( "unimplmented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferGridDirichletOutflowBCVal( const S32 elemIdx, const S32 dim, const VIdx& startVIdx, const VIdx& pdeBufferFaceVSize, REAL& bdryPhi ) {
	/* MODEL START */

	ERROR( "unimplmented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferGridNeumannOutflowBCVal( const S32 elemIdx, const S32 dim, const VIdx& startVIdx, const VIdx& pdeBufferFaceVSize, REAL& bdrySlope ) {
	/* MODEL START */

	ERROR( "unimplmented." );

	/* MODEL END */

	return;
}

void ModelRoutine::adjustPDEBufferGridRHSTimeDependentLinear( const S32 elemIdx, const VIdx& startVIdx, const VIdx& pdeBufferBoxVSize, const REAL phi, REAL& rhs/* INOUT */ ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

void ModelRoutine::updatePDEBufferGridRHSTimeDependentSplitting( const S32 pdeIdx, const VIdx& startVIdx, const VIdx& pdeBufferBoxVSize, const Vector<double>& v_phi/* [idx] */, Vector<double>& v_rhs/* [idx], uptake(-) and secretion (+) */ ) {
	/* MODEL START */

	ERROR( "unimplemented." );

	/* MODEL END */

	return;
}

