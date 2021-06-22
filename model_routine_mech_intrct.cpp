/*

Copyright © 2013, 2017 Battelle Memorial Institute. All Rights Reserved.

NOTICE:  These data were produced by Battelle Memorial Institute (BATTELLE) under Contract No. DE-AC05-76RL01830 with the U.S. Department of Energy (DOE).  For a five year period from May 28, 2013, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, and perform publicly and display publicly, by or on behalf of the Government.  There is provision for the possible extension of the term of this license.  Subsequent to that period or any extension granted, the Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.  The specific term of the license can be identified by inquiry made to BATTELLE or DOE.  NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR BATTELLE, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY DATA, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*/

/* DO NOT USE FUNCTIONS THAT ARE NOT THREAD SAFE (e.g. rand(), use Util::getModelRand() instead) */

#include "biocellion.h"

#include "model_routine.h"

/* MODEL START */

#include "model_define.h"

/* MODEL END */

using namespace std;

#if HAS_SPAGENT
void ModelRoutine::initJunctionSpAgent( const VIdx& vIdx0, const SpAgent& spAgent0, const UBEnv& ubEnv0, const VIdx& vIdx1, const SpAgent& spAgent1, const UBEnv& ubEnv1, const VReal& vDir/* unit direction vector from spAgent1 to spAgent0 */, const REAL dist, BOOL& link, JunctionEnd& end0/* dummy if link == false */, JunctionEnd& end1/* dummy if link == false */ ) {
	/* MODEL START */
    link = false;

    S32 type0 = spAgent0.state.getType();
    S32 type1 = spAgent1.state.getType();
    REAL factor = 1.1;

    REAL R0 = spAgent0.state.getModelReal( CELL_MODEL_REAL_RADIUS ) ;
    REAL R1 = spAgent1.state.getModelReal( CELL_MODEL_REAL_RADIUS );

    REAL dist_threshold = factor * (R0 + R1);

    if ( A_AGENT_BOND_S[type0][type1] > 0.0 ){
 
       if (  dist < dist_threshold )  {
          link = true;
       
          if ( type0 == type1 ) {     
             end0.setType( JUNCTION_END_TYPE_CELL );           
             end1.setType( JUNCTION_END_TYPE_CELL );
          } else {
             end0.setType( JUNCTION_END_TYPE_MICROCARRIER );
             end1.setType( JUNCTION_END_TYPE_MICROCARRIER ); 
          }
       }
       else {
            link = false;
       }
    }


    /* MODEL END */

    return;
}

void ModelRoutine::computeMechIntrctSpAgent( const S32 iter, const VIdx& vIdx0, const SpAgent& spAgent0, const UBEnv& ubEnv0, const VIdx& vIdx1, const SpAgent& spAgent1, const UBEnv& ubEnv1, const VReal& vDir/* unit direction vector from spAgent1 to spAgent0 */, const REAL dist, MechIntrctData& mechIntrctData0, MechIntrctData& mechIntrctData1, BOOL& link, JunctionEnd& end0/* dummy if link == false */, JunctionEnd& end1/* dummy if link == false */, BOOL& unlink  ) {

    /* MODEL START */
    link = false;
    unlink = false;

    S32 type0 = spAgent0.state.getType();
    S32 type1 = spAgent1.state.getType();

    S32 mtype0 = spAgent0.state.getModelInt( CELL_MODEL_INT_STATE ) ;
    S32 mtype1 = spAgent1.state.getModelInt( CELL_MODEL_INT_STATE ) ;

    REAL R0 = A_AGENT_SHOVING_SCALE[mtype0]*spAgent0.state.getModelReal( CELL_MODEL_REAL_RADIUS ) ;
    REAL R1 = A_AGENT_SHOVING_SCALE[mtype1]*spAgent1.state.getModelReal( CELL_MODEL_REAL_RADIUS ) ;
    REAL dist_threshold = R0 + R1;
    REAL D = R0 + R1 - 0.5*A_AGENT_SHOVING_LIMIT[mtype0] - 0.5*A_AGENT_SHOVING_LIMIT[mtype1];
    REAL mag = 0.0;
    REAL stress = 0.0 ; 
    REAL xij  = D - dist ;

    
    if ( A_AGENT_BOND_S[mtype0][mtype1] > 0.0 ){
        REAL sij = A_AGENT_BOND_S[mtype0][mtype1] ;
        if(spAgent0.junctionData.isLinked(spAgent1.junctionData) == true) {
            if( dist > A_AGENT_BOND_DESTROY_FACTOR[mtype0][mtype1]* dist_threshold ) {
                unlink = true;/* break junction */
                //cout << " broken juntion " << endl;
            }
            else{
                // compute elastic force
                REAL Fij = 0.5 * xij * tanh(FABS(xij)*sij);
                mag = mag + Fij ;
                stress = stress + dist * Fij ;

                if ( (type0 == AGENT_CELL_A ) && (type1 == AGENT_MCARRIER ) ){
                    mechIntrctData0.setModelReal( CELL_DIVISION_NORMAL_X,vDir[0] );
                    mechIntrctData0.setModelReal( CELL_DIVISION_NORMAL_Y,vDir[1] );
                    mechIntrctData0.setModelReal( CELL_DIVISION_NORMAL_Z,vDir[2] );
                }
                else if ( (type1 == AGENT_CELL_A ) && (type0 == AGENT_MCARRIER ) ){
                    mechIntrctData1.setModelReal( CELL_DIVISION_NORMAL_X,-vDir[0] );
                    mechIntrctData1.setModelReal( CELL_DIVISION_NORMAL_Y,-vDir[1] );
                    mechIntrctData1.setModelReal( CELL_DIVISION_NORMAL_Z,-vDir[2] );
                } 
 

            }
        }
        else {/* no junction */
            if( dist < A_AGENT_BOND_CREATE_FACTOR[mtype0][mtype1]*dist_threshold ) {

               link = true;/* form junction */

               if ( type0 == type1 ) {
                  end0.setType( JUNCTION_END_TYPE_CELL );
                  end1.setType( JUNCTION_END_TYPE_CELL );
               } else {
                  end0.setType( JUNCTION_END_TYPE_MICROCARRIER );
                  end1.setType( JUNCTION_END_TYPE_MICROCARRIER );
               }

               // add force rigth away
               REAL D = R0 + R1;
               REAL xij  = D - dist  ;
               REAL Fij = 0.5 * xij * tanh(FABS(xij)*sij);
               mag = mag + Fij ;
               stress = stress + dist * Fij;
            }
        }
        
    }
    
    //if ( type0 == type1 ) {
    //   cout << "force " << mag << " xij " << xij << endl;
    //}
    mag = A_AGENT_STIFFNESS[mtype0][mtype1] * mag ;
    stress  = A_AGENT_STIFFNESS[mtype0][mtype1] * stress ;
       

    mechIntrctData0.setModelReal(CELL_MECH_REAL_FORCE_X,vDir[0]*mag);
    mechIntrctData0.setModelReal(CELL_MECH_REAL_FORCE_Y,vDir[1]*mag);
    mechIntrctData0.setModelReal(CELL_MECH_REAL_FORCE_Z,vDir[2]*mag);

    mechIntrctData1.setModelReal(CELL_MECH_REAL_FORCE_X,-vDir[0]*mag);
    mechIntrctData1.setModelReal(CELL_MECH_REAL_FORCE_Y,-vDir[1]*mag);
    mechIntrctData1.setModelReal(CELL_MECH_REAL_FORCE_Z,-vDir[2]*mag); 
 
    mechIntrctData0.setModelReal(CELL_MECH_REAL_STRESS,stress);
    mechIntrctData1.setModelReal(CELL_MECH_REAL_STRESS,stress);


    return;

}
#endif

