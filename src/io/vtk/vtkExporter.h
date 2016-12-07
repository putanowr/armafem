#pragma once

#include <armadillo>
#include <map>
#include <boost/filesystem.hpp>

namespace armafem
{

using FieldsTable = std::map<std::string, arma::vec>;

class VtkExporter
{
  public:
    explicit VtkExporter(const std::string &filename);
    bool exportBar1(const arma::vec &coords, const arma::umat &edof, 
                    const FieldsTable *nodalFields = nullptr,
                    const FieldsTable *cellFields = nullptr);
  private:
    boost::filesystem::path path_;
    std::string title_  = std::string("Generated with use of ArmaFem library");
};

} // namespace armafem
