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
