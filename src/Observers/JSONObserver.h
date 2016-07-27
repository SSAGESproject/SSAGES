#pragma once 

#include "../Simulations/SimObserver.h"
#include "json/json.h"
#include <iostream>
#include <fstream>
#include <memory>

namespace SSAGES
{
	// JSON Observer class that creates a snapshot of the 
	// simulation and dumps it to a JSON file.
	class JSONObserver : public SimObserver
	{
	private:
		std::unique_ptr<std::ofstream> _jsonfs;
		std::string _prefix;
		Json::Value _root;

	public:
		JSONObserver(const std::string& prefix, 
						boost::mpi::communicator world,
						boost::mpi::communicator comm,
						int numwalkers,
						int wid,
						unsigned int frequency = 1);

				// Get Observer name.
		virtual std::string GetName() const override{ return "JSON"; }

		std::string GetPrefix() const {return _prefix;}

		virtual void Visit(const Driver& d) override;

		virtual void Serialize(Json::Value& json) const override;
		
		virtual void PreVisit() override;
		virtual void PostVisit() override;

		~JSONObserver()
		{
			if(_jsonfs)
				_jsonfs->close();
		}
	};
}