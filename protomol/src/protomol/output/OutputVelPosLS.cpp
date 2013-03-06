#include "OutputVelPosLS.h"

using namespace ProtoMol;
const string OutputVelPosLS::keyword("OutputVelPosLS");

OutputVelPosLS::OutputVelPosLS(): OutputVelPos()
{}

OutputVelPosLS::~OutputVelPosLS()
{}

OutputVelPosLS::OutputVelPosLS(const string& filename): OutputVelPos(filename)
{}

