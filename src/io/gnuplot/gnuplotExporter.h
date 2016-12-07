#pragma once

#include <armadillo>
#include <map>
#include <boost/filesystem.hpp>

namespace armafem
{

using FieldsTable = std::map<std::string, arma::vec>;

class GnuplotExporter
{
  public:
    explicit GnuplotExporter(const std::string &path, bool isDataEmbedded=true);
    bool exportBar1(const arma::vec &coords, const arma::umat &edof, 
                    const FieldsTable *nodalFields = nullptr,
                    const FieldsTable *cellFields = nullptr);
    void setTitle(const std::string &title);
    const std::string &getTitle() const;
    void setPath(const std::string &path);
    void setDataEmbedded(bool flag);
  private:
    bool isDataEmbedded_ = true;
    boost::filesystem::path path_;
    std::string title_  = std::string("Generated with use of ArmaFem library");
};

} // namespace armafem
