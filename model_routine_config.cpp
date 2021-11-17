/*

Copyright © 2013, 2017 Battelle Memorial Institute. All Rights Reserved.

NOTICE:  These data were produced by Battelle Memorial Institute (BATTELLE) under Contract No. DE-AC05-76RL01830 with the U.S. Department of Energy (DOE).  For a five year period from May 28, 2013, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, and perform publicly and display publicly, by or on behalf of the Government.  There is provision for the possible extension of the term of this license.  Subsequent to that period or any extension granted, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.  The specific term of the license can be identified by inquiry made to BATTELLE or DOE.  NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR BATTELLE, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY DATA, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*/

/* DO NOT USE FUNCTIONS THAT ARE NOT THREAD SAFE (e.g. rand(), use Util::getModelRand() instead) */

#include "biocellion.h"

#include "model_routine.h"

/* MODEL START */

#include "model_define.h"

#include <string>
// extern "C" {
//   // call to CGAL to create hierarchy of tetrahedra
//   void cfd_setup(char *);
// };

/* MODEL END */

using namespace std;

void ModelRoutine::updateIfGridSpacing( REAL& ifGridSpacing ) {
    /* MODEL START */

    ifGridSpacing = IF_GRID_SPACING;/* sets the grid resolution */

    /* MODEL END */

    return;
}

void ModelRoutine::updateOptModelRoutineCallInfo( OptModelRoutineCallInfo& callInfo ) {
    /* MODEL START */

    callInfo.numComputeMechIntrctIters = 1;
    callInfo.numUpdateIfGridVarPreStateAndGridStepIters = 1;
    callInfo.numUpdateIfGridVarPostStateAndGridStepIters = 0;

    /* MODEL END */

    return;
}

void ModelRoutine::updateDomainBdryType( domain_bdry_type_e a_domainBdryType[DIMENSION] ) {
    /* MODEL START */

    CHECK( DIMENSION == 3 );

    for( S32 dim = 0 ; dim < DIMENSION ; dim++ ) {
        //a_domainBdryType[dim] = DOMAIN_BDRY_TYPE_NONPERIODIC_HARD_WALL;
        a_domainBdryType[dim] = DOMAIN_BDRY_TYPE_PERIODIC;
    }

    /* MODEL END */

    return;
}

void ModelRoutine::updatePDEBufferBdryType( pde_buffer_bdry_type_e& pdeBufferBdryType ) {
    /* MODEL START */

    pdeBufferBdryType = PDE_BUFFER_BDRY_TYPE_HARD_WALL;/* relevant to agents only (not a PDE boundary condition), agents cannot penetrate into the PDE buffer region */

    /* MODEL END */

    return;
}

void ModelRoutine::updateTimeStepInfo( TimeStepInfo& timeStepInfo ) {
    /* MODEL START */

    timeStepInfo.durationBaselineTimeStep = BASELINE_TIME_STEP_DURATION;
    timeStepInfo.numStateAndGridTimeStepsPerBaseline = NUM_STATE_AND_GRID_TIME_STEPS_PER_BASELINE ;

    /* MODEL END */

    return;
}

void ModelRoutine::updateSyncMethod( sync_method_e& mechIntrctSyncMethod, sync_method_e& updateIfGridVarSyncMethod/* dummy if both callUpdateIfGridVarPreStateAndGridStep and callUpdateIfGridVarPostStateAndGridStep are set to false in ModelRoutine::updateOptModelRoutineCallInfo */ ) {
    /* MODEL START */

    mechIntrctSyncMethod = SYNC_METHOD_PER_ATTR;
    updateIfGridVarSyncMethod = SYNC_METHOD_PER_ATTR;

    /* MODEL END */

    return;
}

#if HAS_SPAGENT
void ModelRoutine::updateSpAgentInfo( Vector<SpAgentInfo>& v_spAgentInfo ) {/* set the mechanical interaction range & the numbers of model specific variables */
    /* MODEL START */

    MechModelVarInfo modelVarInfo;

    modelVarInfo.syncMethod = VAR_SYNC_METHOD_DELTA;/* computeMechIntrct */

    v_spAgentInfo.resize( NUM_AGENT_TYPES  );

    for( S32 i = 0 ; i < NUM_AGENT_TYPES; i++ ) {

        SpAgentInfo info;

        info.mechIntrctBdryType = MECH_INTRCT_BDRY_TYPE_SPHERE;
        info.numStateModelReals = NUM_CELL_STATE_REALS;
        info.numStateModelInts = NUM_CELL_STATE_INTS;
        info.numStateInternalModelReals = 0;
        info.numStateInternalModelInts = 0;
        info.v_mechIntrctModelRealInfo.assign( NUM_CELL_MECH_REALS, modelVarInfo );
        info.v_mechIntrctModelIntInfo.clear();
        info.v_boolNetInfo.clear();

        if ( i == AGENT_CELL_A ) {
 
            ODENetInfo odeNetInfo;

            /* ODE setup cell growth */
            odeNetInfo.numVars = NUM_ODE_NET_VAR_GROWING_CELL;
            odeNetInfo.stiff = ODE_STIFF_NORMAL;
            odeNetInfo.h = 0.1;
            odeNetInfo.hm = 0.01;
            odeNetInfo.epsilon = 1e-6;
            odeNetInfo.threshold = 1e-3;
            odeNetInfo.errorThresholdVal = 0.0;
            odeNetInfo.warningThresholdVal = 0.0;
            odeNetInfo.setNegToZero = true;

            info.v_odeNetInfo.push_back( odeNetInfo );
 
        } else {
            info.v_odeNetInfo.clear() ;
        }

        v_spAgentInfo[i] = info;
    }

    /* MODEL END */

    return;
}
#endif

void ModelRoutine::updateJunctionEndInfo( Vector<JunctionEndInfo>& v_junctionEndInfo ) {/* set the numbers of model specific variables */
    /* MODEL START */
    v_junctionEndInfo.resize( NUM_JUNCTION_END_TYPES );
    for( S32 i = 0 ; i < NUM_JUNCTION_END_TYPES ; i++ ) {
        v_junctionEndInfo[i].numModelReals = NUM_JUNCTION_END_REALS;
        v_junctionEndInfo[i].numModelInts = 0;
    }
    /* MODEL END */

    return;
}

void ModelRoutine::updateEllipsoidInfo( EllipsoidInfo& ellipsoidInfo ) {
    /* MODEL START */

    ellipsoidInfo.maxIters = MECH_INTRCT_ELLIPSOID_MAX_ITERS;
    ellipsoidInfo.epsilon = MECH_INTRCT_ELLIPSOID_EPSILON;

    /* MODEL END */

    return;
}

void ModelRoutine::updatePhiPDEInfo( Vector<PDEInfo>& v_phiPDEInfo ) {
    /* MODEL START */

    v_phiPDEInfo.clear();

    /* MODEL END */

    return;
}

void ModelRoutine::updateIfGridModelVarInfo( Vector<IfGridModelVarInfo>& v_ifGridModelRealInfo, Vector<IfGridModelVarInfo>& v_ifGridModelIntInfo, Vector<IfGridModelVarInfo>& v_ifGridInternalModelRealInfo, Vector<IfGridModelVarInfo>& v_ifGridInternalModelIntInfo ) {
    /* MODEL START */
    IfGridModelVarInfo info;
    v_ifGridModelRealInfo.resize(NUM_GRID_MODEL_REALS);
    info.name = "CellDeath";
    info.syncMethod = VAR_SYNC_METHOD_OVERWRITE; // this variable is not expected to be modified by concurrent routines
    v_ifGridModelRealInfo[GRID_MODEL_REAL_DEATHCELL_COUNT ] = info;
    v_ifGridModelIntInfo.clear();
    v_ifGridInternalModelRealInfo.clear();
    v_ifGridInternalModelIntInfo.clear();

    /* MODEL END */

    return;
}

void ModelRoutine::updateRNGInfo( Vector<RNGInfo>& v_rngInfo ) {
    /* MODEL START */

    CHECK( NUM_MODEL_RNGS == 3 );
    v_rngInfo.resize( NUM_MODEL_RNGS );

    RNGInfo rngInfo;

    rngInfo.type = RNG_TYPE_UNIFORM;
    rngInfo.param0 = 0.0;
    rngInfo.param1 = 1.0;
    rngInfo.param2 = 0.0;/* dummy */
    v_rngInfo[ MODEL_RNG_UNIFORM ] =  rngInfo ;

    rngInfo.type = RNG_TYPE_UNIFORM;
    rngInfo.param0 = 0.95;
    rngInfo.param1 = 1.05;
    rngInfo.param2 = 0.0;/* dummy */
    v_rngInfo[MODEL_RNG_UNIFORM_10PERCENT] = rngInfo ;

    rngInfo.type = RNG_TYPE_GAUSSIAN;
    rngInfo.param0 = 0.0;
    rngInfo.param1 = 1.0;
    rngInfo.param2 = 0.0;/* dummy */
    v_rngInfo[MODEL_RNG_GAUSSIAN] = rngInfo;

    /* MODEL END */

    return;
}

void ModelRoutine::updateFileOutputInfo( FileOutputInfo& fileOutputInfo ) {
    /* MODEL START */

    fileOutputInfo.v_spAgentParticleOutput.assign( NUM_AGENT_TYPES , true );
    // CHECK( NUM_PARTICLE_EXTRA_OUTPUT_REALS == 1 );
    fileOutputInfo.v_particleExtraOutputRealName.resize( NUM_PARTICLE_EXTRA_OUTPUT_REALS );
    fileOutputInfo.v_particleExtraOutputRealName[PARTICLE_EXTRA_OUTPUT_REAL_RADIUS] = "radius";
    fileOutputInfo.v_particleExtraOutputRealName[PARTICLE_EXTRA_OUTPUT_REAL_STRESS] = "stress";
    fileOutputInfo.v_particleExtraOutputRealName[PARTICLE_EXTRA_OUTPUT_REAL_ID] = "id";
    fileOutputInfo.v_particleExtraOutputRealName[PARTICEL_EXTRA_OUTPUT_REAL_MicroID] = "microid";
    fileOutputInfo.v_particleExtraOutputRealName[PARTICLE_EXTRA_OUTPUT_REAL_VX] = "vx";
    fileOutputInfo.v_particleExtraOutputRealName[PARTICLE_EXTRA_OUTPUT_REAL_VY] = "vy";
    fileOutputInfo.v_particleExtraOutputRealName[PARTICLE_EXTRA_OUTPUT_REAL_VZ] = "vz";
    fileOutputInfo.v_particleExtraOutputVRealName.clear();
    fileOutputInfo.v_phiOutput.clear();
    fileOutputInfo.v_phiOutputDivideByKappa.clear();

    /* MODEL END */

    return;
}

void ModelRoutine::updateSummaryOutputInfo( Vector<SummaryOutputInfo>& v_summaryOutputRealInfo, Vector<SummaryOutputInfo>& v_summaryOutputIntInfo ) {
    /* MODEL START */

    SummaryOutputInfo info;
    v_summaryOutputIntInfo.clear();
    v_summaryOutputRealInfo.resize( NUM_GRID_SUMMARY_REALS );

    info.name = "Number of microcarriers";
    info.type = SUMMARY_TYPE_SUM;
    v_summaryOutputRealInfo[ GRID_SUMMARY_REAL_MICROCARRIERS  ] = info;

    info.name = "Number of Live Cells";
    info.type = SUMMARY_TYPE_SUM;
    v_summaryOutputRealInfo[GRID_SUMMARY_REAL_LIVE_CELLS] = info;

    info.name = "Number of Live Attached Cells";
    info.type = SUMMARY_TYPE_SUM;
    v_summaryOutputRealInfo[GRID_SUMMARY_REAL_LIVE_ATTACHED] = info;

    info.name = "Number of Death Cells";
    info.type = SUMMARY_TYPE_SUM;
    v_summaryOutputRealInfo[ GRID_SUMMARY_REAL_DEATH ] = info;

    info.name = "Number of Removed Cells";
    info.type = SUMMARY_TYPE_SUM;
    v_summaryOutputRealInfo[ GRID_SUMMARY_REAL_REMOVED ] = info;

    info.name = "Maximum displacement";
    info.type = SUMMARY_TYPE_MAX;
    v_summaryOutputRealInfo[ GRID_SUMMARY_REAL_MAX_DISP ] = info;

    info.name = "Maximum growth rate factor" ;
    info.type = SUMMARY_TYPE_MAX;
    v_summaryOutputRealInfo[ GRID_SUMMARY_REAL_MAX_GROWRATE ] = info;



    /* MODEL END */

    return;
}

void ModelRoutine::initGlobal( Vector<U8>& v_globalData ) {/* called once per simulation */
    /* MODEL START */

    //cfd_setup((char*)"Velocity_220rpm.txt") ;
    
    /* MODEL END */

    return;
}

void ModelRoutine::init( void ) {/* called once per (MPI) process */
    /* MODEL START */

    /* nothing to do */

    /* MODEL END */

    return;
}

void ModelRoutine::term( void ) {/* called once per (MPI) process */
    /* MODEL START */

    /* nothing to do */

    /* MODEL END */

    return;
}

void ModelRoutine::setPDEBuffer( const VIdx& startVIdx, const VIdx& regionVSize, BOOL& isPDEBuffer ) {
    /* MODEL START */

    isPDEBuffer = false;

    /* MODEL END */

    return;
}

void ModelRoutine::setHabitable( const VIdx& vIdx, BOOL& isHabitable ) {
    /* MODEL START */

    isHabitable = true;

    /* MODEL END */

    return;
}

