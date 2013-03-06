#ifndef NORMALMODEOPENMM_H
#define NORMALMODEOPENMM_H

#include <protomol/integrator/openMM/OpenMMIntegrator.h>
#include <protomol/integrator/normal/NormalModeUtilities.h>

namespace ProtoMol {
	class ScalarStructure;
	class ForceGroup;

	class NormalModeOpenMM : public OpenMMIntegrator, public NormalModeUtilities {
		public:
			NormalModeOpenMM();
			NormalModeOpenMM( const std::vector<Value>& base, const std::vector<Value>& params, ForceGroup *forces );
			~NormalModeOpenMM();
			
			virtual std::string getIdNoAlias() const { return keyword; }
			virtual void getParameters( std::vector<Parameter>& parameters ) const;
			virtual unsigned int getParameterSize() const;
			
			virtual void initialize( ProtoMolApp *app );
			virtual void run( int numTimesteps );
		private:
			virtual STSIntegrator *doMake( const std::vector<Value>& values, ForceGroup *fg )const;
		public:
			static const std::string keyword;
		private:
			bool shoudForceRediagOnMinFail;
			bool mRediagOnQuadratic;

			int mRediagonalizationFrequency;
			
			Real mMinimizationLimit;
			int mResiduesPerBlock, mBlockDOF;
			Real mBlockDelta, mSDelta;
			int mModes, mBlockPlatform;
			bool mProtomolDiagonalize;
			NormalModeUtilities *myPreviousNormalMode;
	};
}

#endif
