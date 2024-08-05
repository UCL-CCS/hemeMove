// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <memory>

#include <catch2/catch.hpp>

#include "io/xml.h"
#include "util/Vector3D.h"

#include "tests/helpers/FolderTestFixture.h"

namespace hemelb::tests
{
    namespace xml = hemelb::io::xml;

    namespace {

        xml::Element GetDatum(xml::Element& parent, unsigned iRequired) {
            unsigned i = 0;
            for (xml::Element datum = parent.GetChildOrThrow("datum");
                 datum != xml::Element::Missing();
                 datum = datum.NextSiblingOrNull("datum"))
            {
                if (i == iRequired)
                    return datum;
                ++i;
            }
            throw Exception() << "Cannot find element 'datum' with required index = " << iRequired;
        }

        template <typename T>
        struct TestAttributeTraits {
            static const std::string value;
        };
        template<>
        const std::string TestAttributeTraits<int>::value = "int";
        template<>
        const std::string TestAttributeTraits<unsigned>::value = "unsigned";
        template<>
        const std::string TestAttributeTraits<double>::value = "double";
        template<>
        const std::string TestAttributeTraits<util::Vector3D<int>>::value = "hemelb::util::Vector3D<int>";
        template<>
        const std::string TestAttributeTraits<util::Vector3D<double>>::value = "hemelb::util::Vector3D<double>";

        template<typename T>
        void TestAttributeConversion(int i, const T& EXPECTED, const xml::Element& root) {
            xml::Element shouldWork = root.GetChildOrThrow("conversiontests").GetChildOrThrow("shouldwork");
            xml::Element datum = GetDatum(shouldWork, i);
            auto type = datum.GetAttributeOrThrow("type");
            REQUIRE(TestAttributeTraits<T>::value == type);
            T value;
            datum.GetAttributeOrThrow("value", value);
            REQUIRE(value == EXPECTED);
        }

        template<typename T>
        void TestAttributeConversionFails(int i, const xml::Element& root) {
            xml::Element shouldFail = root.GetChildOrThrow("conversiontests").GetChildOrThrow("shouldfail");
            xml::Element datum = GetDatum(shouldFail, i);
            auto type = datum.GetAttributeOrThrow("type");
            REQUIRE(TestAttributeTraits<T>::value == type);
            T value;
            REQUIRE_THROWS_AS(datum.GetAttributeOrThrow("value", value), io::xml::DeserialisationError);
        }
    }

    TEST_CASE_METHOD(helpers::FolderTestFixture, "XML tests") {
      const std::string testFile{"xmltest.xml"};
      CopyResourceToTempdir(testFile);
      auto xmlDoc = std::make_unique<xml::Document>(testFile);

      SECTION("Test opening a missing file errors") {
          xml::Document xml;
          CHECK_THROWS_AS(xml.LoadFile("this_file_does_not.exist"), xml::ParseError);
      }
      SECTION("TestRead") {
	xml::Element html = xmlDoc->GetRoot();

	// Check that the root elem is html
	REQUIRE(html != xml::Element::Missing());
	REQUIRE(std::string("html") == html.GetName());

	// Check a bunch of elements are present.
	xml::Element head = html.GetChildOrThrow("head");
	{
	  xml::Element title = head.GetChildOrThrow("title");
	}
	xml::Element body = html.GetChildOrThrow("body");
	xml::Element div_banner = body.GetChildOrThrow("div");
	REQUIRE(std::string("banner") == div_banner.GetAttributeOrThrow("id"));
	xml::Element div_header = div_banner.GetChildOrThrow("div");
	REQUIRE(std::string("header") == div_header.GetAttributeOrThrow("id"));

      }
      SECTION("TestSiblings") {
	xml::Element banner = xmlDoc->GetRoot().GetChildOrThrow("body").GetChildOrThrow("div");
	int n = 0;
	std::string expectedIds[] = { "header", "metanav" };

	for (xml::Element div = banner.GetChildOrThrow("div"); div != xml::Element::Missing(); div
	       = div.NextSiblingOrNull("div")) {
	  REQUIRE(expectedIds[n] == div.GetAttributeOrThrow("id"));
	  ++n;
	}
	REQUIRE(2 == n);
      }

      SECTION("TestGetChildNull") {
	xml::Element html = xmlDoc->GetRoot();
	REQUIRE(html.GetChildOrNull("NoSuchElement") == xml::Element::Missing());
      }

      SECTION("TestGetChildThrows") {
	xml::Element html = xmlDoc->GetRoot();
	REQUIRE_THROWS_AS(html.GetChildOrThrow("NoSuchElement"), xml::ChildError);
      }

      SECTION("TestGetParent") {
	xml::Element html = xmlDoc->GetRoot();
	xml::Element body = html.GetChildOrThrow("body");

	xml::Element shouldBeHtml = body.GetParentOrThrow();
	REQUIRE(html == shouldBeHtml);
      }

      SECTION("TestGetParentNull") {
	xml::Element html = xmlDoc->GetRoot();
	REQUIRE(xml::Element::Missing() == html.GetParentOrNull());
      }

      SECTION("TestGetParentThrows") {
	xml::Element html = xmlDoc->GetRoot();
	REQUIRE_THROWS_AS(html.GetParentOrThrow(), xml::ParentError);
      }

        SECTION("Test iteration of named children") {
            auto html = xmlDoc->GetRoot();
            auto n = 0;
            for (auto d: html.Children("body")) {
                n++;
            }
            REQUIRE(n == 1);
        }
        SECTION("Test iteration of all children") {
            auto html = xmlDoc->GetRoot();
            auto n = 0;
            for (auto d: html.Children()) {
                n++;
            }
            REQUIRE(n == 3);
        }

      SECTION("Test attribute conversion") {
	// <datum type="int" value="120" />
	TestAttributeConversion<int>(0, 120, xmlDoc->GetRoot());
	// <datum type="int" value="-24324" />
	TestAttributeConversion<int>(1, -24324, xmlDoc->GetRoot());
	// <datum type="double" value="1.0" />
	TestAttributeConversion<double>(2, 1.0, xmlDoc->GetRoot());
	// <datum type="double" value="1.6e-3" />
	TestAttributeConversion<double>(3, 1.6e-3, xmlDoc->GetRoot());
	// <datum type="hemelb::util::Vector3D<double>" value="(-1.4,11e7,42)" />
	TestAttributeConversion<util::Vector3D<double>>(4, util::Vector3D<double>(-1.4, 11e7, 42), xmlDoc->GetRoot());
	// <datum type="hemelb::util::Vector3D<int>" value="(-1,11,42)" />
	TestAttributeConversion<util::Vector3D<int>>(5, util::Vector3D<int>(-1, 11, 42), xmlDoc->GetRoot());
	// <datum type="unsigned" value="42" />
	TestAttributeConversion<unsigned>(6, 42U, xmlDoc->GetRoot());
      }

      SECTION("Test attributeconversion fails") {
	TestAttributeConversionFails<double>(0, xmlDoc->GetRoot());
	TestAttributeConversionFails<double>(1, xmlDoc->GetRoot());
	TestAttributeConversionFails<int>(2, xmlDoc->GetRoot());
	TestAttributeConversionFails<util::Vector3D<double>>(3, xmlDoc->GetRoot());
	TestAttributeConversionFails<util::Vector3D<double>>(4, xmlDoc->GetRoot());
	TestAttributeConversionFails<util::Vector3D<int>>(5, xmlDoc->GetRoot());
	TestAttributeConversionFails<util::Vector3D<double>>(6, xmlDoc->GetRoot());
	TestAttributeConversionFails<util::Vector3D<int>>(7, xmlDoc->GetRoot());
	TestAttributeConversionFails<unsigned>(8, xmlDoc->GetRoot());
      }

      SECTION("Test that mismatched end tag errors on parse") {
	const std::string testdata = "<sometag>Some text</wrongtag>";
	xml::Document doc;
	REQUIRE_THROWS_AS(doc.LoadString(testdata), xml::ParseError);
      }
    }

}
