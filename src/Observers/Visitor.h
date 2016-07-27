#pragma once

namespace SSAGES
{
	// Forward declare.
	class Driver;

	class Visitor
	{
		// Abstract base class for a visiting object that traverses visitables.
		public:
			virtual void Visit(const Driver& d) = 0;
	};
}
