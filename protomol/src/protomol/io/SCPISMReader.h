#ifndef SCPISMREADER_H
#define SCPISMREADER_H

#include <map>
#include <sstream>
#include <protomol/io/Reader.h>
#include <protomol/topology/CoulombSCPISMParameters.h>

namespace ProtoMol {
    class SCPISMReader : public Reader {
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // Constructors, destructors (both default here), assignment
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        public:
            SCPISMReader();
            explicit SCPISMReader( const std::string &filename );
            virtual ~SCPISMReader();

            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // From class Reader
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        public:
            virtual bool tryFormat();
            virtual bool read();

            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // New methods of class SCPISM
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        public:
            bool read( std::map<std::string, CoulombSCPISMParameters>& inMap );

        private:
            template<typename Type>
            void ConvertString( const std::string& inData, Type& output );
    };
}
#endif // SCPISMREADER_H
