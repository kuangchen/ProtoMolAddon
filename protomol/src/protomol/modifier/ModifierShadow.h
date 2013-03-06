#ifndef MODIFIERSHADOW_H
#define MODIFIERSHADOW_H

#include <protomol/modifier/Modifier.h>
#include <protomol/type/Vector3DBlock.h>

namespace ProtoMol {

    class Integrator;

    //  _________________________________________________________ ModifierShadow

    class ModifierShadow : public Modifier {

        //  ----------------------------------------------------------------  //
        //  Constructors, destructors, assignment
        //  ----------------------------------------------------------------  //
        public:
            ModifierShadow( int order2k, int freq, const Integrator *i,
                            int myId = 10 );


        //  ----------------------------------------------------------------  //
        //  From class Modifier
        //  ----------------------------------------------------------------  //
        private:
            virtual void doExecute(Integrator* i);
            virtual void doInitialize();

            virtual std::string doPrint() const {
                return( std::string( "Modifier Shadow" ) );
            };

        public:
            virtual bool isInternal() const { return( false ); }
    virtual std::string getIdNoAlias() const {return "Shadow";}

        //  ----------------------------------------------------------------  //
        //  New methods of class ModifierShadow.
        //  ----------------------------------------------------------------  //
        public:
            Real calcShadow();
            void fixDiffTable();
            void pauseCalc();
            void resetHistory();   
            void unpauseCalc();

        protected:
            virtual Real getTimestep() const;

        private:
            Real Aij( unsigned int, unsigned int );


        //  ----------------------------------------------------------------  //
        //  My data members
        //  ----------------------------------------------------------------  //
        private:
            const Integrator *myIntegrator;

            Real & myBeta;

            bool paused;

            unsigned int myOrder2k,
                         myShadowK,
                         myFreq;

            std::deque< Vector3DBlock > myStoredPos,
                                        myStoredVel;

            std::deque< Real > myStoredBeta;

            std::vector< Real > coeffs;
            std::vector< std::vector<int> > indeces;


    };

    inline void ModifierShadow::pauseCalc() {
        paused = true;
    }

    inline void ModifierShadow::unpauseCalc() {
        paused = false;
    }


}

#endif /* MODIFIERSHADOW_H */

