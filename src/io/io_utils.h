#pragma once

#include <boost/filesystem.hpp>
#include <string>

namespace armafem {

void fixExtension(boost::filesystem::path pth, const std::string &suffix);

} // namespace armafem
