#include "SimConfig.h"

#include <string>
#include <iostream>
#include <unistd.h>
#include <cstdlib>

#include "util/fileutils.h"
#include "debug/Debugger.h"

namespace hemelb
{

  SimConfig::SimConfig()
  {
    // This constructor only exists to prevent instantiation without
    // using the static load method.
  }

  SimConfig::~SimConfig()
  {
  }

  SimConfig *SimConfig::Load(const char *iPath)
  {
    util::check_file(iPath);
    TiXmlDocument *lConfigFile = new TiXmlDocument();
    lConfigFile->LoadFile(iPath);

    SimConfig *lRet = new SimConfig();

    lRet->DoIO(lConfigFile->FirstChildElement(), true);
    lRet->DataFilePath = util::NormalizePathRelativeToPath(lRet->DataFilePath, iPath);
    delete lConfigFile;

    return lRet;
  }

  void SimConfig::Save(std::string iPath)
  {
    TiXmlDocument lConfigFile;
    TiXmlDeclaration * lDeclaration = new TiXmlDeclaration("1.0", "", "");
    TiXmlElement *lTopElement = new TiXmlElement("hemelbsettings");

    lConfigFile.LinkEndChild(lDeclaration);
    lConfigFile.LinkEndChild(lTopElement);

    DoIO(lTopElement, false);

    lConfigFile.SaveFile(iPath);
  }

  void SimConfig::DoIO(TiXmlElement *iTopNode, bool iIsLoading)
  {
    TiXmlElement* lSimulationElement = GetChild(iTopNode, "simulation", iIsLoading);
    DoIO(lSimulationElement, "cycles", iIsLoading, NumCycles);
    DoIO(lSimulationElement, "cyclesteps", iIsLoading, StepsPerCycle);

    TiXmlElement* lGeometryElement = GetChild(iTopNode, "geometry", iIsLoading);

    if (lGeometryElement != NULL)
    {
      DoIO(GetChild(lGeometryElement, "datafile", iIsLoading), "path", iIsLoading, DataFilePath);
    }

    DoIO(GetChild(iTopNode, "inlets", iIsLoading), iIsLoading, Inlets, "inlet");

    DoIO(GetChild(iTopNode, "outlets", iIsLoading), iIsLoading, Outlets, "outlet");

    TiXmlElement* lVisualisationElement = GetChild(iTopNode, "visualisation", iIsLoading);
    DoIO(GetChild(lVisualisationElement, "centre", iIsLoading), iIsLoading, VisCentre);
    TiXmlElement *lOrientationElement = GetChild(lVisualisationElement, "orientation", iIsLoading);
    DoIO(lOrientationElement, "longitude", iIsLoading, VisLongitude);
    DoIO(lOrientationElement, "latitude", iIsLoading, VisLatitude);

    TiXmlElement *lDisplayElement = GetChild(lVisualisationElement, "display", iIsLoading);

    DoIO(lDisplayElement, "zoom", iIsLoading, VisZoom);
    DoIO(lDisplayElement, "brightness", iIsLoading, VisBrightness);

    TiXmlElement *lRangeElement = GetChild(lVisualisationElement, "range", iIsLoading);

    DoIO(lRangeElement, "maxvelocity", iIsLoading, MaxVelocity);
    DoIO(lRangeElement, "maxstress", iIsLoading, MaxStress);
  }

  void SimConfig::DoIO(TiXmlElement* iParent,
                       std::string iAttributeName,
                       bool iIsLoading,
                       float &value)
  {
    if (iIsLoading)
    {
      char *dummy;
      value = (float) std::strtod(iParent->Attribute(iAttributeName)->c_str(), &dummy);
    }
    else
    {
      // This should be ample.
      char lStringValue[20];

      // %g uses the shorter of decimal / mantissa-exponent notations.
      // 6 significant figures will be written.
      sprintf(lStringValue, "%.6g", value);

      iParent->SetAttribute(iAttributeName, lStringValue);
    }
  }

  void SimConfig::DoIO(TiXmlElement* iParent,
                       std::string iAttributeName,
                       bool iIsLoading,
                       double &value)
  {
    if (iIsLoading)
    {
      const std::string* data = iParent->Attribute(iAttributeName);

      if (data != NULL)
      {
        char *dummy;
        value = std::strtod(iParent->Attribute(iAttributeName)->c_str(), &dummy);
      }
      else
      {
        value = 0.0;
      }
    }
    else
    {
      // This should be ample.
      char lStringValue[20];

      // %g uses the shorter of decimal / mantissa-exponent notations.
      // 6 significant figures will be written.
      sprintf(lStringValue, "%.6g", value);

      iParent->SetAttribute(iAttributeName, lStringValue);
    }
  }

  void SimConfig::DoIO(TiXmlElement* iParent,
                       std::string iAttributeName,
                       bool iIsLoading,
                       std::string &iValue)
  {
    if (iIsLoading)
    {
      if (iParent->Attribute(iAttributeName) == 0)
        iValue = "";
      else
        iValue = std::string(iParent->Attribute(iAttributeName)->c_str());
    }
    else
    {
      if (iValue != "")
        iParent->SetAttribute(iAttributeName, iValue);
    }
  }

  void SimConfig::DoIO(TiXmlElement* iParent,
                       std::string iAttributeName,
                       bool iIsLoading,
                       long &bValue)
  {
    if (iIsLoading)
    {
      char *dummy;
      // Read in, in base 10.
      bValue = std::strtol(iParent->Attribute(iAttributeName)->c_str(), &dummy, 10);
    }
    else
    {
      // This should be ample.
      char lStringValue[20];

      // %ld specifies long integer style.
      sprintf(lStringValue, "%ld", bValue);

      iParent->SetAttribute(iAttributeName, lStringValue);
    }
  }

  void SimConfig::DoIO(TiXmlElement* iParent,
                       std::string iAttributeName,
                       bool iIsLoading,
                       unsigned long &bValue)
  {
    if (iIsLoading)
    {
      char *dummy;
      // Read in, in base 10.
      bValue = std::strtoul(iParent->Attribute(iAttributeName)->c_str(), &dummy, 10);
    }
    else
    {
      // This should be ample.
      char lStringValue[20];

      // %ld specifies long integer style.
      sprintf(lStringValue, "%ld", bValue);

      iParent->SetAttribute(iAttributeName, lStringValue);
    }
  }

  void SimConfig::DoIO(TiXmlElement *iParent,
                       bool iIsLoading,
                       std::vector<InOutLet> &bResult,
                       std::string iChildNodeName)
  {
    if (iIsLoading)
    {
      TiXmlElement *lCurrentLet = iParent->FirstChildElement(iChildNodeName);

      while (lCurrentLet != NULL)
      {
        InOutLet lNew;
        DoIO(lCurrentLet, iIsLoading, lNew);
        bResult.push_back(lNew);
        lCurrentLet = lCurrentLet->NextSiblingElement(iChildNodeName);
      }
    }
    else
    {
      for (unsigned int ii = 0; ii < bResult.size(); ii++)
      {
        // NB we're good up to 99 io-lets here.
        DoIO(GetChild(iParent, iChildNodeName, iIsLoading), iIsLoading, bResult[ii]);
      }
    }
  }

  void SimConfig::DoIO(TiXmlElement *iParent, bool iIsLoading, InOutLet &value)
  {
    TiXmlElement* lPositionElement = GetChild(iParent, "position", iIsLoading);
    TiXmlElement* lNormalElement = GetChild(iParent, "normal", iIsLoading);
    TiXmlElement* lPressureElement = GetChild(iParent, "pressure", iIsLoading);

    DoIO(lPressureElement, "path", iIsLoading, value.PFilePath);

    // TODO In the case of a sinusoidal pressure condition, the specification of minimum and maximum
    // is redundant and therefore should be removed. In the case of a file, they should be
    // calculated from the file anyway.

    if (value.PFilePath == "")
    {
      DoIO(lPressureElement, "mean", iIsLoading, value.PMean);
      DoIO(lPressureElement, "amplitude", iIsLoading, value.PAmp);
      DoIO(lPressureElement, "phase", iIsLoading, value.PPhase);
    }
    DoIO(lPressureElement, "minimum", iIsLoading, value.PMin);
    DoIO(lPressureElement, "maximum", iIsLoading, value.PMax);

    DoIO(lPositionElement, iIsLoading, value.Position);
    DoIO(lNormalElement, iIsLoading, value.Normal);
  }

  void SimConfig::DoIO(TiXmlElement *iParent, bool iIsLoading, Vector &iValue)
  {
    DoIO(iParent, "x", iIsLoading, iValue.x);
    DoIO(iParent, "y", iIsLoading, iValue.y);
    DoIO(iParent, "z", iIsLoading, iValue.z);
  }

  TiXmlElement *SimConfig::GetChild(TiXmlElement *iParent,
                                    std::string iChildNodeName,
                                    bool iIsLoading)
  {
    if (iIsLoading)
    {
      return iParent->FirstChildElement(iChildNodeName);
    }
    else
    {
      TiXmlElement* lNewChild = new TiXmlElement(iChildNodeName);
      iParent->LinkEndChild(lNewChild);
      return lNewChild;
    }
  }

}
