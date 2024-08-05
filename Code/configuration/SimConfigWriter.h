// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_CONFIGURATION_SIMCONFIGWRITER_H
#define HEMELB_CONFIGURATION_SIMCONFIGWRITER_H

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "configuration/SimConfig.h"
#include "io/xml.h"

namespace hemelb::extraction { struct PropertyOutputFile; }

namespace hemelb::configuration
{
    struct MonitoringConfig;

    // Turn a SimConfig into an XML file
    class SimConfigWriter {
        using path = std::filesystem::path;
        using Element = io::xml::Element;
        using Document = io::xml::Document;

        path outputXmlPath;
        std::unique_ptr<Document> outputXml;

    public:
        explicit SimConfigWriter(path p);

        virtual ~SimConfigWriter() = default;

        virtual void Write(const SimConfig&);

        // Turn a full path into and XML-relative path
        [[nodiscard]] std::string FullPathToRelPath(const path& p) const;

        virtual void DoIOForSimulation(const GlobalSimInfo&);
        virtual void DoIOForGeometry(const path&);

        virtual void DoIOForInOutlets(std::string type, std::vector<IoletConfig> const&);

        void DoIOForBaseInOutlet(Element ioletEl, IoletConfigBase const& ioletConf) const;

        void DoIOForCosinePressureInOutlet(Element&, CosinePressureIoletConfig const&) const;
        void DoIOForFilePressureInOutlet(Element&, FilePressureIoletConfig const&) const;
        void DoIOForMultiscalePressureInOutlet(Element&, MultiscalePressureIoletConfig const&) const;
        void DoIOForParabolicVelocityInOutlet(Element&, ParabolicVelocityIoletConfig const&) const;
        void DoIOForWomersleyVelocityInOutlet(Element&, WomersleyVelocityIoletConfig const&) const;
        void DoIOForFileVelocityInOutlet(Element&, FileVelocityIoletConfig const&) const;

        virtual void DoIOForProperties(std::vector<extraction::PropertyOutputFile> const &) const;
        virtual void DoIOForInitialConditions(ICConfig const&);
        virtual void DoIOForMonitoring(const MonitoringConfig& mon_conf);

//        virtual Element TemplateCellConfig readCell(const Element& cellNode) const;
//        virtual Element std::map<std::string, TemplateCellConfig> readTemplateCells(Element const& cellsEl) const;
//        virtual Element RBCConfig DoIOForRedBloodCells(SimConfig const&, const Element& rbcEl) const;

    private:
        template <typename... CStrs>
        bool InOriginal(CStrs... path) const;

        void AddPathChild(Element&, const path&) const;

    };

}

#endif
