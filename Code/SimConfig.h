#ifndef HEMELB_SIMCONFIG_H
#define HEMELB_SIMCONFIG_H

#include <string.h>
#include <vector>
#include "xml/tinyxml.h"
#include "util/Vector3D.h"
#include "lb/boundaries/InOutLets.h"

namespace hemelb
{
  class SimConfig
  {
    public:
      static SimConfig *Load(const char *iPath);
      ~SimConfig();

      void Save(std::string iPath);

      std::string DataFilePath;
      std::vector<lb::boundaries::InOutLet*> Inlets;
      std::vector<lb::boundaries::InOutLet*> Outlets;
      util::Vector3D VisCentre;
      float VisLongitude;
      float VisLatitude;
      float VisZoom;
      float VisBrightness;
      float MaxVelocity;
      float MaxStress;
      unsigned long NumCycles;
      long StepsPerCycle;

      void DoIO(TiXmlElement *iXmlNode, bool iIsLoading, lb::boundaries::InOutLetCosine* value);
      void DoIO(TiXmlElement *iXmlNode, bool iIsLoading, lb::boundaries::InOutLetFile* value);
    private:
      SimConfig();
      void DoIO(TiXmlElement *iXmlNode, bool iIsLoading);
      void DoIO(TiXmlElement* iXmlNode, std::string iAttributeName, bool iIsLoading, long &value);
      void DoIO(TiXmlElement* iXmlNode,
                std::string iAttributeName,
                bool iIsLoading,
                unsigned long &value);
      void DoIO(TiXmlElement* iXmlNode, std::string iAttributeName, bool iIsLoading, float &value);
      void DoIO(TiXmlElement* iXmlNode, std::string iAttributeName, bool iIsLoading, double &value);
      void DoIO(TiXmlElement *iXmlNode,
                std::string iAttributeName,
                bool iIsLoading,
                std::string &iValue);
      void DoIO(TiXmlElement *iXmlNode,
                bool iIsLoading,
                std::vector<lb::boundaries::InOutLet*> &value,
                std::string iChildNodeName);
      void DoIO(TiXmlElement *iXmlNode, bool iIsLoading, util::Vector3D &iValue);
      TiXmlElement* GetChild(TiXmlElement *iParent, std::string iChildNodeName, bool iIsLoading);
  };
}

#endif /* HEMELB_SIMCONFIG_H */
