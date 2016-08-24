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

namespace SSAGES
{
	// Forward declare.
	class Driver;

    //! Base class for an object that traverses visitables.
    /*!
     * Abstract base class for a visiting object that traverses visitables.
     *
     * \ingroup Core
     */
	class Visitor
	{
		public:
            //! Visit
            /*!
             * \param d Driver to be visited.
             */
			virtual void Visit(const Driver& d) = 0;

            //! Destructor.
            virtual ~Visitor() {}
	};
}
