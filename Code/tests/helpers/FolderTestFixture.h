// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_HELPERS_FOLDERTESTFIXTURE_H
#define HEMELB_TESTS_HELPERS_FOLDERTESTFIXTURE_H

#include <string>
#include <filesystem>
#include <fstream>

#include "tests/helpers/HasCommsTestFixture.h"

namespace hemelb::io::xml { class Document; }

namespace hemelb::tests::helpers
{
    void ModifyXMLInput(io::xml::Document &document, std::vector<std::string> const& elements,
                        std::string const& _value);

    //! \brief Modify XML document
    //! \details HemeLB parameters cannot be modified programmatically, so we have to jump
    //! these hoops to test it.
    //! \param[in] document: Document to modify
    //! \param[in] elements: hierarchy of elements + attribute (last item)
    //!   Should not include "hemelbsettings"
    //! \param[in] value: Value to set the attribute to
    template<class T>
    void ModifyXMLInput(io::xml::Document &document, std::vector<std::string> const& elements,
                        T const &_value)
    {
        std::ostringstream attr_value;
        attr_value << _value;
        ModifyXMLInput(document, elements, attr_value.str());
    }

    //! \brief Modify XML document by deleting an element if it exists
    //! \details HemeLB parameters cannot be modified programmatically, so we have to jump
    //! these hoops to test it.
    //! \param[in] document: Document to modify
    //! \param[in] elements: hierarchy of elements, last item will be removed.
    //!   Should not include "hemelbsettings"
    //! \param[in] value: Value to set the attribute to
    void DeleteXMLInput(io::xml::Document &document, std::vector<std::string> const& elements);

    class FolderTestFixture : public HasCommsTestFixture
    {
    private:
        using path = std::filesystem::path;
        path origin;
        path tempPath;
        bool ownsTempPath = false;

    public:
        FolderTestFixture();
        ~FolderTestFixture();

    protected:
        void ReturnToOrigin();
        void CopyResourceToTempdir(const std::string & resource) const;
        //! \brief Modify XML resource
        //! \details HemeLB parameters cannot be modified programmatically, so we have to jump
        //! these hoops to test it.
        //! \param[in] resource: Name of the resource.
        //!   Should exist in the temp dir prior to call.
        //! \param[in] elements: hierarchy of elements + attribute (last item)
        //!   Should not include "hemelbsettings"
        //! \param[in] value: Value to set the attribute to
        void ModifyXMLInput(std::string const &resource, std::vector<std::string> const& elements,
                            std::string const &_value) const;
        template<class T>
        void ModifyXMLInput(std::string const &resource, std::vector<std::string> const& elements,
                            T const &_value) const
        {
            std::ostringstream attr_value;
            attr_value << _value;
            ModifyXMLInput(resource, elements, attr_value.str());
        }

        [[nodiscard]] path ConstructTempPath() const;
        void DeleteXMLInput(std::string const &resource, std::vector<std::string> const& elements) const;
        void MoveToTempdir();
        void AssertPresent(const std::string &fname);
        const path& GetTempdir();
    };
}
#endif
