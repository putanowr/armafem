#include "io/io_utils.h"

namespace bf = boost::filesystem;

namespace armafem {

void fixExtension(boost::filesystem::path pth, const std::string &suffix)
{
   if (bf::extension(pth) != suffix)
   {
     pth += suffix;
   }
}

} // namespace armafem
