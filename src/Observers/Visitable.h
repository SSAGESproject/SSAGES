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
