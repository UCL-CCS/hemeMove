// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REDBLOOD_MESHIO_H
#define HEMELB_REDBLOOD_MESHIO_H

#include <memory>
#include <string>
#include <vector>

#include <redblood/Mesh.h>

class vtkPolyData;
template <class> class vtkSmartPointer;

namespace hemelb::util {
    class UnitConverter;
}

namespace hemelb::redblood {
    class CellBase;
    class Cell;

    // Base class for reading/writing meshes.
    // 
    // Use non-virtual interface pattern, so concrete derived classes
    // only have to implement two member functions (read and write).
    // 
    class MeshIO {
    public:
      using MeshPtr = std::shared_ptr<MeshData>;
      using ConstMeshPtr = std::shared_ptr<MeshData>;
      using ScalarField = std::pair<std::string, std::vector<double>>;
      using PointScalarData = std::vector<ScalarField>;

      // For UnitConverter arguments, use nullptr (i.e. absent) to
      // indicate I/O in lattice units.

      // Where are we writing to?
      enum class Storage {file, string};

      virtual ~MeshIO() = default;

      // Read mesh from file or string.
      [[nodiscard]] MeshPtr readFile(std::string const &filename, bool fixFacetOrientation) const;
      [[nodiscard]] MeshPtr readString(std::string const &data, bool fixFacetOrientation) const;

      // Overload set to write mesh to a file in physical units.
      // First arg: file name
      // Second arg: mesh as either `MeshData` or (`Vertices`, `Facets`)
      // Third/Fourth arg: UnitConverter
      // Fourth/Fifth arg: optional scalar data to attach to the points
      void writeFile(std::string const &filename,
		     MeshData const& mesh,
		     util::UnitConverter const*,
		     PointScalarData const& data = {}) const;
      void writeFile(std::string const &filename,
		     MeshData::Vertices const &, MeshData::Facets const &,
		     util::UnitConverter const*,
		     PointScalarData const& data = {}) const;
      // Overload set to write mesh to a string in physical units (return value).
      // Mesh as above
      // Unit converter
      // Optional scalar data
      [[nodiscard]] std::string writeString(MeshData const& mesh,
			      util::UnitConverter const*,
			      PointScalarData const& data = {}) const;
      [[nodiscard]] std::string writeString(MeshData::Vertices const &, MeshData::Facets const &,
			      util::UnitConverter const*,
			      PointScalarData const& data = {}) const;

      //! Write cell-mesh to file or string
      void writeFile(std::string const& filename, CellBase const&, util::UnitConverter const*) const;
      [[nodiscard]] std::string writeString(CellBase const&, util::UnitConverter const*) const;

      //! Write cell-mesh including individual forces to file or string
      void writeFile(std::string const& filename, Cell const&, util::UnitConverter const*) const;
      std::string writeString(Cell const&, util::UnitConverter const*) const;

    private:
      // Derived class implement these to actually do the serialisation.
      [[nodiscard]] virtual MeshPtr read(Storage,
                           std::string const& filename_or_data,
                           bool fixFacetOrientation) const = 0;
      // Discard is fine iff storage is file
      virtual std::string write(Storage m, std::string const& filename,
                                MeshData::Vertices const &, MeshData::Facets const &,
                                util::UnitConverter const* c, PointScalarData const& data) const = 0;
    };

    // Format is from T. Krueger's thesis
    class KruegerMeshIO : public MeshIO {
    public:
      ~KruegerMeshIO() override = default;

      [[nodiscard]] MeshPtr read(Storage,
                   std::string const &filename_or_data,
                   bool fixFacetOrientation) const override;

      std::string write(Storage m,
                        std::string const &filename,
                        MeshData::Vertices const &, MeshData::Facets const &,
                        util::UnitConverter const* c,
                        PointScalarData const& data) const override;
    };

    // VTK XML PolyData format
    class VTKMeshIO : public MeshIO {
    public:
      using PolyDataPtr = vtkSmartPointer<vtkPolyData>;
      // Helper member function, exposed to allow testing of the fix orientations
      [[nodiscard]] std::tuple<MeshPtr, PolyDataPtr> readUnoriented(Storage, std::string const &) const;

      ~VTKMeshIO() override = default;

      [[nodiscard]] MeshPtr read(Storage,
                   std::string const &filename_or_data,
                   bool fixFacetOrientation) const override;

      std::string write(Storage m,
                        std::string const &filename,
                        MeshData::Vertices const &, MeshData::Facets const &,
                        util::UnitConverter const* c,
                        PointScalarData const& data) const override;

    };

}
#endif
