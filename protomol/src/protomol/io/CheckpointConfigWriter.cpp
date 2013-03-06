#include <protomol/io/CheckpointConfigWriter.h>

using namespace ProtoMol;

CheckpointConfigWriter::CheckpointConfigWriter() : Writer() {}

CheckpointConfigWriter::CheckpointConfigWriter(const std::string &filename) :
  Writer(filename) {}

bool CheckpointConfigWriter::write(const int &id, const int &steps,
                                   const Random &rand,
                                   const Integrator *integ) {
  file
    << "!Checkpoint File!" << std::endl
    << "#ID" << std::endl
    << id << std::endl
    << "#Step" << std::endl
    << steps << std::endl
    << "#Random" << std::endl
    << rand << std::endl
    << "#Integrator" << std::endl
    << *integ << std::endl;

  return !file.fail();
}
