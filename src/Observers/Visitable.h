/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Ben Sikora <bsikora906@gmail.com>
 *
 * SSAGES is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SSAGES is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SSAGES.  If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

#include "Visitor.h"

namespace SSAGES
{
	//! Base class for visitable objects.
    /*
     * Abstract base class for visitable objects, traversed (usually) by loggers.
     *
     * \ingroup Core
     */
	class Visitable
	{
		public:
            //! Accept a visitor
            /*
             * \param v Reference to the object visiting.
             *
             * This function is called when a visitor (usually a looger) wants
             * to visit this object.
             */
			virtual void AcceptVisitor(Visitor &v) const = 0;

            //! Destructor.
            virtual ~Visitable() {}
	};
}
