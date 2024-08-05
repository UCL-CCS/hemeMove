// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_REDBLOOD_RBCINSERTER_H
#define HEMELB_REDBLOOD_RBCINSERTER_H

#include <iostream>
#include <memory>
#include <functional>
#include <list>
#include <random>
#include <utility>

#include "io/xml.h"
#include "lb/iolets/InOutLet.h"
#include "log/Logger.h"
#include "redblood/types.h"
#include "redblood/Mesh.h"
#include "redblood/Cell.h"
#include "units.h"
#include "util/Matrix3D.h"

namespace hemelb::redblood
{

    /**
     * The RBC Inserter inserts cells into the simulation while some condition
     * evaluates to true.  The cell is inserted into an inlet or the start of a
     * flow extension (if the inlet has a flow extension defined).
     * The shape of the cell is read from a text file.
     */
    class RBCInserter
    {
    public:
        /**
         * Creates an RBC Inserter.
         *
         * @param condition a cell will only be inserted on a LB step if this
         * condition evaluates to true
         * @param cell to insert. Inserted as is. E.g. same position etc.
         * @param scale the scale of the cell to insert
         */
        RBCInserter(std::function<bool()> condition, std::unique_ptr<CellBase const> cell) :
            condition(std::move(condition)), cell(std::move(cell)), barycentre(this->cell->GetBarycentre())
        {
        }
        RBCInserter(RBCInserter &&c) :
            condition(std::move(c.condition)), cell(std::move(c.cell)),
                barycentre(c.barycentre)
        {
        }
        RBCInserter(RBCInserter const&c) :
            condition(c.condition), cell(c.cell->clone()), barycentre(c.barycentre)
        {
        }
        virtual ~RBCInserter() = default;

        void operator=(RBCInserter const &c)
        {
          condition = c.condition;
          cell = c.cell->clone();
        }
        void operator=(RBCInserter &&c)
        {
          condition = std::move(c.condition);
          cell = std::move(c.cell);
        }

        /**
         * Cell insertion callback called on each step of the simulation.  Cells
         * are inserted into the inlets while condition() evaluates to true.
         * Multiple cells are potentially inserted on each iteration.
         *
         * @param insertFn the function to insert a new cell into the simulation
         *
         * @see hemelb::redblood::CellArmy::SetCellInsertion
         * @see hemelb::redblood::CellArmy::CallCellInsertion
         */
        void operator()(CellInserter insertFn)
        {
          if (condition())
          {
            log::Logger::Log<log::Debug, log::OnePerCore>("Dropping one cell at (%f, %f, %f)",
                                                          barycentre.x(),
                                                          barycentre.y(),
                                                          barycentre.z());
            insertFn(drop());
          }
        }

        virtual CellContainer::value_type drop()
        {
          return CellContainer::value_type(cell->clone().release());
        }

      private:
        //! When to insert cells
        std::function<bool()> condition;
        //! The shape of the cells to insert
        std::unique_ptr<CellBase const> cell;
        //! barycentre -- for logging
        LatticePosition barycentre;
    };

    //! Red blood cell inserter that adds random rotation and translation to each cell
    class RBCInserterWithPerturbation : public RBCInserter
    {
    public:
        using PRNG = std::minstd_rand;

        RBCInserterWithPerturbation(std::function<bool()> condition,
                                    std::unique_ptr<CellBase const> cell,
                                    util::Matrix3D const &initialRotation, Angle dtheta, Angle dphi,
                                    LatticePosition const& dx, LatticePosition const& dy,
                                    PRNG::result_type randomGeneratorSeed = PRNG::default_seed) :
            RBCInserter(condition, std::move(cell)), initialRotation(initialRotation),
                dtheta(dtheta), dphi(dphi), dx(dx), dy(dy), randomGenerator(randomGeneratorSeed), uniformDistribution(-1.0,1.0)
        {
        }

        //! Rotates and translates the template cell according to random dist
        CellContainer::value_type drop() override;

        using RBCInserter::operator();
      private:
        //! Rotation to flow axis + offset
        util::Matrix3D initialRotation;
        //! Range of values over which to pick extra rotation.
        Angle dtheta, dphi;
        //! Two vectors alongst which to move cell
        LatticePosition dx, dy;

        PRNG randomGenerator;
        std::uniform_real_distribution<double> uniformDistribution;
    };

    //! Composite RBCInserter that inserts multiple cells at each LB step
    class CompositeRBCInserter
    {
      public:

        void AddInserter(const std::shared_ptr<RBCInserter> & inserter)
        {
          inserters.push_back(inserter);
        }

        void RemoveInserter(const std::shared_ptr<RBCInserter> & inserter)
        {
          auto pos = std::find(std::begin(inserters), std::end(inserters), inserter);
          if (pos != std::end(inserters))
          {
            inserters.erase(pos);
          }
        }

        std::size_t Size() const noexcept
        {
          return inserters.size();
        }

        void operator()(CellInserter insertFn) const
        {
          for (const std::shared_ptr<RBCInserter> & inserter : inserters)
          {
            (*inserter)(insertFn);
          }
        }
      private:
        std::list<std::shared_ptr<RBCInserter>> inserters;
    };
}

#endif
