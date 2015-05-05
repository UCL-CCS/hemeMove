#include "io/xml/XmlAbstractionLayer.h"
#include "RBCInserter.h"

namespace hemelb
{
  namespace redblood
  {

RBCInserter::RBCInserter(std::function<bool()> condition,
                         const std::string & mesh_path,
                         std::vector<lb::iolets::InOutLet *> inlets,
                         Cell::Moduli moduli,
                         Dimensionless scale) :
    condition(condition), inlets(inlets), moduli(moduli), scale(scale)
{
  this->shape = read_mesh(mesh_path);
}

RBCInserter::RBCInserter(std::function<bool()> condition,
                         std::istream & mesh_stream,
                         std::vector<lb::iolets::InOutLet *> inlets,
                         Cell::Moduli moduli,
                         Dimensionless scale) :
    condition(condition), inlets(inlets), moduli(moduli), scale(scale)
{
  this->shape = read_mesh(mesh_stream);
}

RBCInserter::RBCInserter(std::function<bool()> condition,
                         const MeshData & shape,
                         std::vector<lb::iolets::InOutLet *> inlets,
                         Cell::Moduli moduli,
                         Dimensionless scale) :
    condition(condition), shape(new MeshData(shape)), inlets(inlets), moduli(moduli), scale(scale)
{
}

void RBCInserter::operator()(std::function<void(CellContainer::value_type)> insertFn) const
{
  for (auto inlet: inlets) {
    if (condition()) {
      std::shared_ptr<Cell> cell = std::make_shared<Cell>(this->shape->vertices,
                                                          Mesh(*this->shape),
                                                          this->scale);
      cell->moduli = this->moduli;
      const std::shared_ptr<FlowExtension> flowExt = inlet->GetFlowExtension();
      *cell += inlet->GetPosition();
      if (flowExt)
        *cell += flowExt->normal * flowExt->length;
      insertFn(cell);
    }
    else
      break;
  }
}

void RBCInserter::SetShape(const MeshData & shape)
{
  this->shape.reset(new MeshData(shape));
}
std::shared_ptr<const MeshData> RBCInserter::GetShape() const
{
  return this->shape;
}

void RBCInserter::AddInLet(lb::iolets::InOutLet * inlet)
{
  this->inlets.push_back(inlet);
}
void RBCInserter::RemoveInLet(lb::iolets::InOutLet * inlet)
{
  this->inlets.erase(
      std::remove(std::begin(this->inlets), std::end(this->inlets), inlet),
      std::end(this->inlets));
}

void RBCInserter::SetScale(Dimensionless scale)
{
  this->scale = scale;
}
Dimensionless RBCInserter::GetScale() const
{
  return this->scale;
}

void RBCInserter::SetModuli(Cell::Moduli & moduli)
{
  this->moduli = moduli;
}
const Cell::Moduli & RBCInserter::GetModuli() const {
  return this->moduli;
}

void RBCInserter::SetCondition(std::function<bool()> condition)
{
  this->condition = condition;
}
std::function<bool()> RBCInserter::GetCondition() const
{
  return this->condition;
}

  }
}
