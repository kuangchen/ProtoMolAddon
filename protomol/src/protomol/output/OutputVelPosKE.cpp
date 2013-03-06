#include "OutputVelPosKE.h"

using namespace ProtoMol;
const string OutputVelPosKE::keyword("OutputVelPosKE");

OutputVelPosKE::OutputVelPosKE(): OutputVelPos()
{}

OutputVelPosKE::~OutputVelPosKE()
{}

OutputVelPosKE::OutputVelPosKE(const string& filename): OutputVelPos(filename)
{}

