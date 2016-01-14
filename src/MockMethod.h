#pragma once 

#include "Method.h"
#include <iostream>

namespace SSAGES
{
	class MockMethod : public Method
	{
	public:
		MockMethod(unsigned int frequency) : 
		Method(frequency)
		{
		}

		void PreSimulation(Snapshot* snapshot, const CVList& cvs) override
		{
		}

		void PostIntegration(Snapshot* snapshot, const CVList& cvs) override
		{
		}

		void PostSimulation(Snapshot* snapshot, const CVList& cvs) override
		{
		}
	};
}