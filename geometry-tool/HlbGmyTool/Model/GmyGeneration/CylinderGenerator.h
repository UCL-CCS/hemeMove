// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HLBGMYTOOL_GMY_CYLINDERGENERATOR_H
#define HLBGMYTOOL_GMY_CYLINDERGENERATOR_H

#include "GeometryGenerator.h"

#include "GenerationError.h"
#include "Index.h"

namespace hemelb::gmytool::gmy {

class GeometryWriter;
class Site;
class BlockWriter;

struct CylinderData {
  // Cylinder parameters
  Vector Centre;
  Vector Axis;
  double Radius;
  double Length;
};

class CylinderGenerator : public GeometryGenerator {
 public:
  CylinderGenerator();
  virtual ~CylinderGenerator();

  inline void SetCylinderCentre(Vector v) { this->Cylinder->Centre = v; }
  inline void SetCylinderAxis(Vector n) { this->Cylinder->Axis = n; }
  inline void SetCylinderRadius(double r) { this->Cylinder->Radius = r; }
  inline void SetCylinderLength(double l) { this->Cylinder->Length = l; }

 private:
  virtual void ComputeBounds(double[]) const;
  void ClassifySite(Site& site);
  void ComputeCylinderNormalAtAPoint(Vector& wallNormal,
                                     const Vector& surfacePoint,
                                     const Vector& cylinderAxis) const;
  CylinderData* Cylinder;

 protected:
  virtual int BlockInsideOrOutsideSurface(const Block& block) { return 0; }
};

}  // namespace hemelb::gmytool::gmy

#endif  // HLBGMYTOOL_GMY_CYLINDERGENERATOR_H
