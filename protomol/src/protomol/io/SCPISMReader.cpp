#include <protomol/io/SCPISMReader.h>

#include <protomol/base/StringUtilities.h>
#include <protomol/base/Report.h>
#include <protomol/base/MathUtilities.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

SCPISMReader::SCPISMReader() :
        Reader() {

}

SCPISMReader::SCPISMReader( const string &filename ) :
        Reader( filename ) {
}

SCPISMReader::~SCPISMReader() {

}

bool SCPISMReader::tryFormat() {
    if ( !open() ) {
        return false;
    }

    return file.good();
}

bool SCPISMReader::read(){
    return !file.fail();
}

bool SCPISMReader::read( std::map<std::string, CoulombSCPISMParameters>& inMap ) {
    if ( !tryFormat() ) {
        return false;
    }

    if ( !open() ) {
        return false;
    }

    while ( !file.eof() ) {
        std::string line;
        std::getline( file, line );

        int endPos = line.find( '!' );
        if ( endPos != -1 ) {
            line = line.substr( 0, endPos );
        }

        if ( line.substr( 0, 3 ) != "end" ) {
            CoulombSCPISMParameters tempPerameter;

            std::string atomName;
            ConvertString( line.substr( 0, 4 ), atomName );

            double alpha_i = 0.0;
            ConvertString( line.substr( 5, 8 ), alpha_i );

            double hbond_factor = 0.0;
            ConvertString( line.substr( 13, 8 ), hbond_factor );

            double R_iw = 0.0;
            ConvertString( line.substr( 21, 8 ), R_iw );

            double unused = 0.0;
            ConvertString( line.substr( 29, 8 ), unused );

            double r_cov = 0.0;
            ConvertString( line.substr( 37, 8 ), r_cov );

            double gamma_i = 0.0;
            ConvertString( line.substr( 45, 7 ), gamma_i );

            std::string bond;
            ConvertString( line.substr( 52, 3 ), bond );

            Hbonded eBond;

            if ( bond == "PH" ){
                eBond = PH;
            }else if ( bond == "PA" ){
                eBond = PA;
            }else{
                eBond = NO;
            }

            tempPerameter.set(
                atomName,
                alpha_i,
                hbond_factor,
                R_iw,
                r_cov,
                gamma_i,
                eBond
            );

            inMap[ atomName ] = tempPerameter;
        } else {
            break;
        }
    }

    return !file.fail();
}

template<typename Type>
void SCPISMReader::ConvertString( const std::string& inData, Type& output ) {
    std::istringstream stream( inData );

    stream >> output;
}
