#ifndef IOInterface_hpp__
#define IOInterface_hpp__

#include <fstream>
#include <vector>
#include <math.h>
#include <iostream>
#include <map>
#include <cstring>

namespace evo
{
    class IOInterface
    {
    public:

        void set_fname(std::string file);                                                     /* setting IO file name */
        void fgetdata( size_t samples, size_t variants, std::vector<std::vector<int>> &out ); /* Reads a binary SNP data stored in the PLINK bed format */
        void fgetdata( bool includes_id, std::vector<std::vector<int>> &out );                /* Reads ASCII SNP data with extra info (the very first two columns -> outputs in the additional data structure) */
        void fgetdata( std::vector<std::vector<int>> &out );                                  /* Reads ASCII SNP data */
        template <typename T>
        void fgetdata( std::vector<std::vector<T>> &out );                                /* Reads general ASCII formated data */

    private:

        void str_parse( std::string& snpStr, std::vector<int>& markers );                     /* Parsing a SNP string */
        std::string io_file;                                                                  /* IO file name */
        std::map<size_t, size_t> snp_id;                                                      /* KEY: the consecutive index; the VALUE: observation ID */
    };

} // end of namespace evo

#endif // IOInterface_hpp__
