// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "tests/helpers/FolderTestFixture.h"

#include <sstream>
#include <fstream>
#include <iomanip>
#include <filesystem>

#include <unistd.h>

#include <catch2/catch.hpp>

#include "io/xml.h"
#include "resources/Resource.h"
#include "util/clock.h"

namespace fs = std::filesystem;

namespace hemelb::tests::helpers
{
    using Document = io::xml::Document;
    using Element = io::xml::Element;

    void ModifyXMLInput(Document &document, std::vector<std::string> const& elements,
                        std::string const& _value)
    {
        auto child = document.GetRoot();
        if (child.GetName() != "hemelbsettings")
            throw (Exception() << "Root element not called 'hemelbsettings'");


        // Point to the actual last *element* (cos last name in vector is the attribute).
        std::string const& attribute = elements.back();
        auto end = --elements.end();
        // Use iterators to avoid last val
        for (auto iter = elements.begin(); iter != end; ++iter) {
            auto& name = *iter;
            auto next_child = child.GetChildOrNull(name.c_str());
            if (next_child  == Element::Missing()) {
                next_child = child.AddChild(name.c_str());
            }
            child = next_child;
        }
        child.SetAttribute(attribute.c_str(), _value.c_str());
    }

    //! \brief Modify XML document by deleting an element if it exists
    //! \details HemeLB parameters cannot be modified programmatically, so we have to jump
    //! these hoops to test it.
    //! \param[in] document: Document to modify
    //! \param[in] elements: hierarchy of elements, last item will be removed.
    //!   Should not include "hemelbsettings"
    //! \param[in] value: Value to set the attribute to
    void DeleteXMLInput(Document &document, std::vector<std::string> const& elements)
    {
        auto element = document.GetRoot();
        if (element.GetName() != "hemelbsettings")
            throw (Exception() << "Root element not called 'hemelbsettings'");
        for (auto& name : elements) {
            auto next_child = element.GetChildOrNull(name.c_str());
            if (next_child  == Element::Missing())
                return;

            element = next_child;
        }
        element.GetParentOrThrow().DeleteChild(element);
    }

    FolderTestFixture::FolderTestFixture()
    {
        auto comm = Comms();
        // We will do the MPI comms using std::string instead of std::filesystem::path
        std::string temp_string;
        if (comm.OnIORank()) {
            tempPath = ConstructTempPath();
            // create a folder to work in
            fs::create_directory(tempPath);
            ownsTempPath = true;
            temp_string = tempPath;
        }
        comm.Broadcast(temp_string, comm.GetIORank());
        if (!comm.OnIORank()) {
            tempPath = temp_string;
        }

        // store current location
        origin = fs::current_path();
        MoveToTempdir();
    }

    FolderTestFixture::~FolderTestFixture() {
        ReturnToOrigin();
        if (ownsTempPath) {
            fs::remove_all(tempPath);
        }
    }

    namespace {
        std::int64_t HackyUUID() {
            static std::int32_t uuid[2] = {getpid(), std::int32_t(util::clock() * 100000)};
            ++uuid[1];
            return *reinterpret_cast<std::int64_t*>(uuid);
        }
    }

    auto FolderTestFixture::ConstructTempPath() const -> path {
        std::stringstream tmp_name;
        tmp_name << "HemeLBTest-" << std::fixed << HackyUUID();
        // TODO: find a portable uuid solution. BOOST?
        return fs::temp_directory_path() / std::move(tmp_name).str();
    }

    void FolderTestFixture::ReturnToOrigin()
    {
        // return to origin
        fs::current_path(origin);
    }

    void FolderTestFixture::CopyResourceToTempdir(const std::string & resource) const
    {
        bool ok = fs::copy_file(resources::Resource(resource).Path(),
                                tempPath / resource);
        REQUIRE(ok);
    }

    void FolderTestFixture::ModifyXMLInput(std::string const &resource, std::vector<std::string> const& elements,
                                           std::string const &_value) const
    {
        auto const filename = tempPath / resource;
        Document document;
        document.LoadFile(filename);
        helpers::ModifyXMLInput(document, elements, _value);
        document.SaveFile(filename);
    }

    void FolderTestFixture::DeleteXMLInput(std::string const &resource, std::vector<std::string> const& elements) const
    {
        std::string const filename = tempPath / resource;
        Document document;
        document.LoadFile(filename);
        helpers::DeleteXMLInput(document, elements);
        document.SaveFile(filename);
    }

    void FolderTestFixture::MoveToTempdir()
    {
        fs::current_path(tempPath);
    }

    void FolderTestFixture::AssertPresent(const std::string &fname)
    {
        // "does directory exist" actually works for files too.
        REQUIRE(fs::exists(fname));
    }
    const fs::path& FolderTestFixture::GetTempdir()
    {
        return tempPath;
    }
}

