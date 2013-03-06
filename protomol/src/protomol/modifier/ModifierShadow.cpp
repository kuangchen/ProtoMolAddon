#include <protomol/modifier/ModifierShadow.h>
#include <protomol/integrator/Integrator.h>
#include <protomol/topology/Topology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/base/PMConstants.h>
#include <protomol/ProtoMolApp.h>

using namespace ProtoMol::Report;
using std::deque;
using std::vector;

namespace ProtoMol {

    //  _________________________________________________________ ModifierShadow

    ModifierShadow::ModifierShadow( int order2k, int freq,
                                    const Integrator *intg, int myId ) :
        Modifier(myId), myIntegrator(intg), myBeta( myIntegrator->myBeta ),
        paused(false), myOrder2k(order2k), myShadowK(0), myFreq(freq) {
                          
        //  ----------------------------------------------------------------  //
        //  Check for valid order of shadow Hamiltonian.  2k must equal one   //
        //  of these values.                                                  //
        //  ----------------------------------------------------------------  //

        switch( myOrder2k ) {

            case 4:
            case 8:
            case 12:
            case 16:
            case 20:
            case 24:
                myShadowK = myOrder2k >> 1;
                break;

            default:
                report << error
                       << "Incorrect value for parameter shadowEnergy: ("
                       << myOrder2k << "). "
                       << "Currently only implemented for multiples of 4."
                       << endr;
                break;

        }


    }


    //  --------------------------------------------------------------------  //
    //  --------------------------------------------------------------------  //

    void ModifierShadow::doInitialize() {

        //  ----------------------------------------------------------------  //
        //  Initialize coeffs with for A_{ij}'s.  Taken from Engle et al.,    //
        //  2005, "Monitoring Energy Drift With Shadow Hamiltonians".         //
        //  ----------------------------------------------------------------  //

        switch( myShadowK ) {

            case 2: {
                Real cfs[3] = { 1., -1./2, 2./3 };
                coeffs = vector< Real >( cfs, cfs + 3 ); }
                break;

            case 4: {
                Real cfs[10] = { 1., -3./2, 16./7, 13./21, -32./21, 36./35,
                                 -5./84, 22./105, -9./35, 4./35 };
                coeffs = vector< Real >( cfs, cfs + 10 ); }
                break;

            case 6: {
                Real cfs[21] = { 1., -5./2, 54./11, 74./33, -72./11, 125./22,
                                 -19./22, 141./44, -375./88, 500./231, 29./220,
                                 -3./5, 325./308, -200./231, 45./154, -7./1320,
                                 137./4620, -125./1848, 5./63, -15./308, 1./77 };

                coeffs = vector< Real >( cfs, cfs + 21 ); }
                break;

            case 8: {
                Real cfs[36] = { 1., -7./2, 128./15, 73./15, -256./15,
                                 10976./585, -41./12, 1712./117, -2744./117,
                                 10976./715, 743./585, -3712./585,
                                 406112./32175, -43904./3575, 6860./1287,
                                 -31./130, 4896./3575, -104272./32175,
                                 128968./32175, -3430./1287, 3136./3861,
                                 37./1925, -28544./225225, 11368./32175,
                                 -1568./2925, 140./297, -896./3861, 112./2145,
                                 -761./1801800, 22./6825, -343./32175,
                                 1918./96525, -175./7722, 28./1755, -14./2145,
                                 8./6435 };

                coeffs = vector< Real >( cfs, cfs + 36 ); }
                break;

            case 10: {
                Real cfs[55] = { 1., -9./2, 250./19, 484./57, -2000./57,
                                 30375./646, -497./57, 167125./3876,
                                 -212625./2584, 21600./323, 34167./6460,
                                 -9625./323, 88695./1292, -25920./323,
                                 185220./4199, -14917./7752, 46975./3876,
                                 -83025./2584, 192600./4199, -154350./4199,
                                 666792./46189, 2761./6783, -19300./6783,
                                 501525./58786, -417600./29393, 652680./46189,
                                 -381024./46189, 110250./46189, -831./18088,
                                 4925./13832, -564975./470288, 43740./19019,
                                 -6615./2431, 186543./92378, -165375./184756,
                                 9000./46189, 4861./2116296, -10475./529074,
                                 14985./198968, -53520./323323, 21315./92378,
                                 -882./4199, 875./7106, -2000./46189,
                                 675./92378, -671./21162960, 7129./23279256,
                                 -6849./5173168, 99./29393, -1029./184756,
                                 2877./461890, -875./184756, 10./4199,
                                 -135./184756, 5./46189 };
                    
                coeffs = vector< Real >( cfs, cfs + 55 ); }
                break;

            case 12: {
                Real cfs[78] = { 1., -11./2, 432./23, 905./69, -1440./23,
                                 15972./161, -1635./92, 702./7, -35937./161,
                                 665500./3059, 12126./805, -76896./805,
                                 3865224./15295, -1064800./3059,
                                 24257475./104006, -953./115, 126376./2185,
                                 -378004./2185, 6355525./22287, -8085825./29716,
                                 6899904./52003, 45454./15295, -345888./15295,
                                 19489833./260015, -7320500./52003,
                                 118053045./728042, -41399424./364021,
                                 4024944./96577, -8315./12236, 1168245./208012,
                                 -8509083./416024, 62523725./1456084,
                                 -331518825./5824336, 230715540./4732273,
                                 -2515590./96577, 705672./96577,
                                 176005./1872108, -131564./156009,
                                 3674891./1092063, -25688300./3276189,
                                 222509925./18929092, -55582560./4732273,
                                 752136./96577, -313632./96577, 136125./193154,
                                 -2339./328440, 126393./1820105,
                                 -1106061./3640210, 1311035./1670214,
                                 -50132115./37858184, 2117016./1391845,
                                 -30492./25415, 307098./482885, -81675./386308,
                                 24200./676039, 58301./240253860,
                                 -51684./20021155, 588181./47322730,
                                 -506990./14196819, 98901./1456084,
                                 -2119392./23661365, 40194./482885,
                                 -26136./482885, 2475./104006, -4400./676039,
                                 594./676039, (-6617./6597360)/437,
                                 (83711./7147140)/437, (-81191./1299480)/437,
                                 78419./170361828, -75339./75716368,
                                 35937./23661365, -1617./965770, 4521./3380195,
                                 -4125./5408312, 605./2028117, -99./1352078,
                                 6./676039 };

                coeffs = vector< Real >( cfs, cfs + 78 ); }
                break;

            default:
                report << error
                       << "Incorrect k parameter, " << myShadowK
                       << ", for shadow Hamiltonian.  "
                       << "Currently only implemented for even 'k'."
                       << endr;
                break;

        }


        //  ----------------------------------------------------------------  //
        //  Initialize A_{ij} indeces based on predefined coefficients.       //
        //  Follows pattern A_10, A_20, A_21, A_30, A_31, A_32, A_40 ...      //
        //  These correspond to the ith and jth backwards differences.        //
        //  ----------------------------------------------------------------  //

        unsigned int count = 0,
                     i = 0;

        while( count < coeffs.size() ) {

            i++;

            for( unsigned int j = 0; j < i; j++ ) {

                count++;

                vector< int > index;

                index.push_back(i);
                index.push_back(j);

                indeces.push_back( index );

            }

        }


        //  ----------------------------------------------------------------  //
        //  ----------------------------------------------------------------  //

        resetHistory();

    }


    //  --------------------------------------------------------------------  //
    //  --------------------------------------------------------------------  //

    void ModifierShadow::fixDiffTable() {

        //  ----------------------------------------------------------------  //
        //  ----------------------------------------------------------------  //

        Vector3DBlock tempBlock( app->positions.size() );

        unsigned int size = myStoredBeta.size() - 1;


        for( unsigned int i = 0; i < size; i++ ) {

            for( unsigned int j = size; j > i; j-- ) {

                myStoredBeta[j] = myStoredBeta[j-1] - myStoredBeta[j];

                tempBlock.intoAssign( myStoredPos[j-1] );
                tempBlock.intoSubtract( myStoredPos[j] );
                myStoredPos[j].intoAssign( tempBlock );

                tempBlock.intoAssign( myStoredVel[j-1] );
                tempBlock.intoSubtract( myStoredVel[j] );
                myStoredVel[j].intoAssign( tempBlock );

            }

        }

    }


    //  --------------------------------------------------------------------  //
    //  --------------------------------------------------------------------  //

    Real ModifierShadow::getTimestep() const {

        return( myIntegrator->getTimestep() );

    }


    //  --------------------------------------------------------------------  //
    //  --------------------------------------------------------------------  //

    void ModifierShadow::doExecute(Integrator* i) {

        //  FIXME: frequency does not give correct answer.
        //  if( myFreq ) ;

        if( paused )
            return;

        app->energies[ScalarStructure::SHADOW] = calcShadow();

    }


    //  --------------------------------------------------------------------  //
    //  Caculate the Aij.                                                     //
    //  --------------------------------------------------------------------  //

    Real ModifierShadow::Aij( unsigned int i, unsigned int j ) {

        Real a_ij = 0.;

        //  ----------------------------------------------------------------  //
        //  ----------------------------------------------------------------  //

        for( unsigned int k = 0; k < app->positions.size(); k++ ) {

            a_ij += (myStoredPos[i])[k].dot( (myStoredVel[j])[k] ) *
                    app->topology->atoms[k].scaledMass;

            a_ij -= (myStoredVel[i])[k].dot( (myStoredPos[j])[k] ) *
                    app->topology->atoms[k].scaledMass;

        }


        //  ----------------------------------------------------------------  //
        //  ----------------------------------------------------------------  //

        //Beta x Alpha, Alpha = 0 if j\ne 0, 1 otherwise;
        //  FIXME
        if( j == 0 )
            a_ij -= myStoredBeta[i];

        if( i == 0 )
            a_ij += myStoredBeta[j];

        a_ij /= ( 2. * getTimestep() * Constant::INV_TIMEFACTOR );

        return( a_ij );

    }


    //  --------------------------------------------------------------------  //
    //  --------------------------------------------------------------------  //

    void ModifierShadow::resetHistory() {

        //  Reset the queues and beta term.                                     
        myStoredBeta.clear();
        myStoredPos.clear();
        myStoredVel.clear();

        myBeta = 0.;

    }


    //  --------------------------------------------------------------------  //
    //  Calculate the Shadow Hamiltonian based on Bob Skeel/David Hardy's     //
    //  interpolated method, using Dan Engle's backward differences.          //
    //  --------------------------------------------------------------------  //

    Real ModifierShadow::calcShadow() {

        //  ----------------------------------------------------------------  //
        //  Store current values.  If running reverse timesteps, we can       // 
        //  ignore this function for now.                                     //
        //  ----------------------------------------------------------------  //

        if( myIntegrator->isForward() ) {

            myStoredBeta.push_front( myBeta );
            myStoredPos.push_front( app->positions );
            myStoredVel.push_front( app->velocities );

        }
        else
            return( 0. );


        //  ----------------------------------------------------------------  //
        //  Since we add the current values to the queue, we will eventually  //
        //  exceed the number we need to store.                               //
        //  ----------------------------------------------------------------  //

        while( myStoredBeta.size() > ( myShadowK + 1 ) ) {

            myStoredBeta.pop_back();
            myStoredPos.pop_back();
            myStoredVel.pop_back();

        }


        //  ----------------------------------------------------------------  //
        //  Recompute the backwards differences with the current values       // 
        //  stored in position [0].  Formula goes: value[i] = value[i-1] -    //
        //  value[i], where value[0] are the current values. Formula is       //
        //  defined in EnSD05, section 5.1.                                   // 
        //  ----------------------------------------------------------------  //

        Vector3DBlock tempBlock( app->positions.size() );

        for( unsigned int i = 1; i < myStoredBeta.size(); i++ ) {

            myStoredBeta[i] = myStoredBeta[i-1] - myStoredBeta[i];

            tempBlock.intoAssign( myStoredPos[i-1] );
            tempBlock.intoSubtract( myStoredPos[i] );
            myStoredPos[i].intoAssign( tempBlock );

            tempBlock.intoAssign( myStoredVel[i-1] );
            tempBlock.intoSubtract( myStoredVel[i] );
            myStoredVel[i].intoAssign( tempBlock );

        }


        //  ----------------------------------------------------------------  //
        //  If we do not yet have enough stored values to correctly compute   //
        //  the shadow, return. We need k + 1 values.                         //
        //  ----------------------------------------------------------------  //

        if( myStoredBeta.size() < myShadowK + 1 )
            return( 0. );


        //  ----------------------------------------------------------------  //
        //  Calculate the shadow.                                             //
        //  ----------------------------------------------------------------  //

        Real shadow = 0.;

        for( unsigned int i = 0; i < coeffs.size(); i++ )
            shadow += coeffs[i] * Aij( indeces[i][0], indeces[i][1] );


        //  TODO: Remove after testing.
        /*
        Real total = myEnergies->potentialEnergy() +
                     kineticEnergy( myTopology, myVelocities );

        fprintf( stderr, "dbg: %10s = %22.18lf", "shadow", shadow );
        fprintf( stderr, "%10s = %22.18lf", "total", total );
        fprintf( stderr, "%10s = %22.18lf\n", "diff", shadow - total );
        */

        return( shadow );

    }


    //  --------------------------------------------------------------------  //

}

