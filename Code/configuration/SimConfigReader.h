// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_CONFIGURATION_SIMCONFIGREADER_H
#define HEMELB_CONFIGURATION_SIMCONFIGREADER_H

#include <filesystem>
#include <map>
#include <string>
#include <string_view>
#include <vector>

#include "Exception.h"
#include "quantity.h"
#include "configuration/MonitoringConfig.h"
#include "configuration/SimConfig.h"
#include "extraction/OutputField.h"
#include "extraction/PropertyOutputFile.h"
#include "io/xml.h"

namespace hemelb::extraction {
    class PlaneGeometrySelector;
    class StraightLineGeometrySelector;
    class SurfacePointSelector;
}

namespace hemelb::configuration {

    void CheckNoAttributes(const io::xml::Element& elem);
    void CheckNoChildren(const io::xml::Element& elem);
    void CheckEmpty(const io::xml::Element& elem);

    // Note on specifying physical quantities.
    //
    // In the XML, these are encoded in attributes of an element (the name
    // of the element doesn't matter for these purposes). The two required
    // attributes are:
    // 1. "units" (e.g. "m" for a physical length or "lattice" for something specified in scaled units)
    // 2. "value" which is a string encoding of the representation (typically a double)
    //
    // Now some parts of the code (e.g. RBC) can specify units as multiple options (e.g. physical OR lattice)
    // for these, use the quantity union types. The consumer (i.e. SimBuilder etc) will have
    // to take care around this.
    inline void check_unit_spec(const io::xml::Element& elem, std::string_view actual, std::string_view const& expected) {
        if (actual != expected)
            throw Exception() << "Invalid units for element " << elem.GetPath()
                              << "."" Expected '" << expected
                              << "', got '" << actual << "'";
    }

    //! Check the units of the quantity and decode the value into @param value
    template<typename T>
    void PopDimensionalValue(io::xml::Element& elem, std::string_view units, T& value)
    {
        auto got = elem.PopAttributeOrThrow("units");
        check_unit_spec(elem, got, units);
        elem.PopAttributeOrThrow("value", value);
        CheckEmpty(elem);
    }

    //! Check the units of the quantity and return the value
    template<typename T>
    T PopDimensionalValue(io::xml::Element& elem, std::string_view units)
    {
        auto got = elem.PopAttributeOrThrow("units");
        check_unit_spec(elem, got, units);
        auto ans = elem.PopAttributeOrThrow<T>("value");
        CheckEmpty(elem);
        return ans;
    }

    //! Given an element (@param elem), check for a child with the given @param name.
    //! If it exists, return the unit-checked (against @param unit) value.
    //! If it doesn't exist, return @param default_value.
    template <typename T>
    T PopDimensionalValueWithDefault(io::xml::Element& el,
                                     char const* name, char const* unit, T default_value) {
        if (el)
            if (auto qEl = el.PopChildOrNull(name))
                return PopDimensionalValue<T>(*qEl, unit);
        return default_value;
    }

    template <QuantityUnion QUnion, IsVariantAlternative<QUnion> DefaultQ>
    QUnion PopDimensionalValueWithDefault(io::xml::Element& elem, char const* name, DefaultQ const& default_q) {
        using RepT = typename quantity_union_traits<QUnion>::representation_type;
        using io::xml::Element;
        if (elem)
            if (auto qEl = elem.PopChildOrNull(name)) {
                auto xml_units = qEl->PopAttributeOrThrow("units");
                auto xml_val = qEl->PopAttributeOrThrow<RepT>("value");
                CheckEmpty(*qEl);
                return quantity_union_factory<QUnion>()(xml_val, xml_units);
            }
        return default_q;
    }

    class SimConfigReader {
        using path = std::filesystem::path;
        using Element = io::xml::Element;
        path xmlFilePath;

    public:
        SimConfigReader(path);

        virtual ~SimConfigReader() = default;

        [[nodiscard]] virtual SimConfig Read() const;

        // Turn an input XML-relative path into a full path
        [[nodiscard]] path RelPathToFullPath(std::string_view path) const;

        /**
         * Check that the iolet is OK for the CMake configuration.
         * @param ioletEl
         * @param requiredBC
         */
        virtual void CheckIoletMatchesCMake(const Element &ioletEl,
                                            std::string_view requiredBC) const;

        [[nodiscard]] virtual SimConfig DoIO(Element xmlNode) const;

        [[nodiscard]] virtual GlobalSimInfo DoIOForSimulation(Element simEl) const;

        [[nodiscard]] virtual path DoIOForGeometry(Element geometryEl) const;

        [[nodiscard]] virtual std::vector <IoletConfig>
        DoIOForInOutlets(GlobalSimInfo const &sim_info, Element xmlNode) const;

        void
        DoIOForBaseInOutlet(GlobalSimInfo const &sim_info, Element &ioletEl, IoletConfigBase &ioletConf) const;

        [[nodiscard]] IoletConfig DoIOForPressureInOutlet(Element &ioletEl) const;

        [[nodiscard]] IoletConfig DoIOForCosinePressureInOutlet(Element &conditionEl) const;

        [[nodiscard]] IoletConfig DoIOForFilePressureInOutlet(Element &conditionEl) const;

        [[nodiscard]] IoletConfig DoIOForMultiscalePressureInOutlet(
                Element &conditionEl) const;

        [[nodiscard]] IoletConfig DoIOForVelocityInOutlet(Element &ioletEl) const;

        [[nodiscard]] IoletConfig DoIOForParabolicVelocityInOutlet(
                Element &conditionEl) const;

        /**
         * Reads a Womersley velocity iolet definition from the XML config file and returns
         * an InOutLetWomersleyVelocity object
         *
         * @param ioletEl in memory representation of <inlet> or <outlet> xml element
         * @return InOutLetWomersleyVelocity object
         */
        [[nodiscard]] IoletConfig DoIOForWomersleyVelocityInOutlet(Element &conditionEl) const;

        /**
         * Reads a file velocity iolet definition from the XML config file and returns
         * an InOutLetFileVelocity object
         *
         * @param ioletEl in memory representation of <inlet> or <outlet> xml element
         * @return InOutLetFileVelocity object
         */
        [[nodiscard]] IoletConfig DoIOForFileVelocityInOutlet(Element &conditionEl) const;

        [[nodiscard]] virtual std::vector <extraction::PropertyOutputFile>
        DoIOForProperties(GlobalSimInfo const &sim_info, const Element &xmlNode) const;

        [[nodiscard]] virtual extraction::OutputField
        DoIOForPropertyField(GlobalSimInfo const &sim_info, Element &xmlNode) const;

        [[nodiscard]] virtual extraction::PropertyOutputFile
        DoIOForPropertyOutputFile(GlobalSimInfo const &sim_info, Element &propertyoutputEl) const;

        [[nodiscard]] extraction::StraightLineGeometrySelector *DoIOForLineGeometry(
                Element &xmlNode) const;

        [[nodiscard]] extraction::PlaneGeometrySelector *DoIOForPlaneGeometry(Element &) const;

        [[nodiscard]] extraction::SurfacePointSelector *DoIOForSurfacePoint(Element &) const;

        [[nodiscard]] virtual ICConfig DoIOForInitialConditions(Element parent) const;
        //virtual void DoIOForCheckpointFile(const Element& checkpointEl) const;

        /**
         * Reads monitoring configuration from XML file
         *
         * @param monEl in memory representation of <monitoring> xml element
         */
        [[nodiscard]] virtual MonitoringConfig DoIOForMonitoring(Element &monEl) const;

        /**
         * Reads configuration of steady state flow convergence check from XML file
         *
         * @param convEl in memory representation of the <steady_flow_convergence> XML element
         */
        void DoIOForSteadyFlowConvergence(Element &convEl, MonitoringConfig &monitoringConfig) const;

        /**
         * Reads the configuration of one of the possible several converge criteria provided
         *
         * @param criterionEl in memory representation of the <criterion> XML element
         */
        void DoIOForConvergenceCriterion(Element &criterionEl, MonitoringConfig &monitoringConfig) const;

        [[nodiscard]] virtual TemplateCellConfig readCell(Element &cellNode) const;

        [[nodiscard]] virtual std::map <std::string, TemplateCellConfig>
        readTemplateCells(Element& cellsEl) const;

        [[nodiscard]] virtual RBCConfig DoIOForRedBloodCells(SimConfig const &, Element &rbcEl) const;

    };
}
#endif