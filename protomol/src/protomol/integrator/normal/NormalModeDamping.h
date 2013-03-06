/*  -*- c++ -*-  */
#ifndef NORMALMODEDAMPING_H
#define NORMALMODEDAMPING_H

#include <protomol/integrator/STSIntegrator.h>
#include <protomol/integrator/normal/NormalModeUtilities.h>

#include <protomol/io/DCDTrajectoryReader.h>
#include <protomol/io/XYZTrajectoryWriter.h>


namespace ProtoMol {

    class ScalarStructure;
    class ForceGroup;

    //__________________________________________________ NormalModeDamping
    class NormalModeDamping : public STSIntegrator, public NormalModeUtilities {
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // Constructors, destructors, assignment
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    public:
        NormalModeDamping();
        NormalModeDamping(Real timestep, int firstmode, int nummode,
            std::string dcdf, std::string modeff, std::string atomff, Real temp,
            ForceGroup *overloadedForces);
        ~NormalModeDamping(); 

        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // New methods of class NormalModeDamping
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    protected:
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // From class Makeable
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    public:
        virtual std::string getIdNoAlias() const{return keyword;}
        virtual unsigned int getParameterSize() const{return 7;}
        virtual void getParameters(std::vector<Parameter>& parameters) const;

        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // From class Integrator
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    public:
        virtual void initialize(ProtoMolApp* appp);
        virtual void run(int numTimesteps);
    protected:
        //virtual void addModifierAfterInitialize();

        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // From class STSIntegrator
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    private:
        virtual STSIntegrator* doMake(const std::vector<Value>& values, ForceGroup* fg)const;

        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // New methods of class NormalModeDamping
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    public:
        //virtual void forceProjection();
    private:	
        void modeProjector(Vector3DBlock *inp, double *tmpc, bool forcep);

        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // My data members
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    public:
        static const std::string keyword;

    private:
        Real randStp, myTemperature;
        int aveForceCount, counter;
        // dcd input and x0
        std::string myDCDFile;
        // mode force output:
        std::string modeForceFile, atomForceFile; 

        DCDTrajectoryReader myReaderF, myReaderX, myReaderV;
        Vector3DBlock *waterForces, x0;

        XYZTrajectoryWriter *myWriter;

        double *posC, *velC;
        
        ofstream myFile;

    };
}

#endif


