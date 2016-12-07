#include "io/gnuplot/gnuplotExporter.h"
#include "io/io_utils.h"

using namespace arma;

namespace armafem
{

GnuplotExporter::GnuplotExporter(const std::string &path, bool isDataEmbedded) : path_(path), isDataEmbedded_(isDataEmbedded)
{
  fixExtension(path_, "gpt");
}

bool GnuplotExporter::exportBar1(const vec &coords, const umat &edof, 
                             const FieldsTable *nodalFields,
                             const FieldsTable *cellFields)
{
  std::ofstream out;
  out.open(path_.string().c_str());
  if (nullptr != nodalFields)
  {
    auto nf = nodalFields->size();
    out << "plot \\" << std::endl;
    for (auto item : *nodalFields)
    {
      out << "'-' w lp title(\"" << item.first <<"\"), \\" << std::endl;
    }
    out << std::endl;
    auto ndof = coords.size();
    auto field = mat(ndof, 2, fill::zeros);
    field(span::all, 0) = coords;
    for (auto item : *nodalFields)
    {
      field(span::all, 1) = item.second;
      out << field;
      out << "e" << std::endl;
    }
  }
}

} // namespace armafem
