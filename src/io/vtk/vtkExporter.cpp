#include "io/vtk/vtkExporter.h"
#include "io/io_utils.h"

using namespace arma;

namespace armafem
{

VtkExporter::VtkExporter(const std::string &filename) : path_(filename)
{
  fixExtension(path_, "vtk");
}

bool VtkExporter::exportBar1(const vec &coords, const umat &edof, 
                             const FieldsTable *nodalFields,
                             const FieldsTable *cellFields)
{
  std::ofstream out;
  out.open(path_.string().c_str());
  out << "# vtk DataFile Version 2.0" << std::endl;
  out << title_ << std::endl;
  out << "ASCII" << std::endl;
  out << "DATASET UNSTRUCTURED_GRID" << std::endl;

  out << "POINTS " << coords.size() << " double" << std::endl;
  for (auto x : coords) out << x << " 0.0 0.0" <<  std::endl;
  auto cells = edof;
  auto nc = edof.n_rows;
  cells(span::all,0).fill(2);
  out << "CELLS " << nc << " " << nc*3 << std::endl;
  out << cells << std::endl;
  out << "CELL_TYPES " << nc << std::endl;
  for (decltype(nc) i=0; i<nc; ++i) out << 3 << std::endl;
}

} // namespace armafem
