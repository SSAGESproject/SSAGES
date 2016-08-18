#include "JSONObserver.h"
#include "../Drivers/Driver.h"
#include <ctime>
#include <iomanip>
#include <stdio.h>

namespace SSAGES
{
	//! Local function, undocumented
	/*!
	 * \return Timestamp
	 */
	std::string GetTimestamp()
	{
		time_t rawtime;
		struct tm * timeinfo;
		char buffer[80];

		time (&rawtime);
		timeinfo = localtime(&rawtime);

  		strftime(buffer,80,"%Y-%m-%d %r",timeinfo);
  		return std::string(buffer);
	}

	JSONObserver::JSONObserver(const std::string& prefix,
								boost::mpi::communicator world,
								boost::mpi::communicator comm,
								int numwalkers,
								int wid,
								unsigned int frequency) : 
	SimObserver(world, comm, numwalkers, wid, frequency), _prefix(prefix), _root(), _writetoother(true)
	{
		if(_comm.rank() != 0)
			return;

		std::string Name = _prefix + "_" + std::to_string(_wid) +".json";
		std::string BckupName = Name + "_backup";
		
		_jsonfs = std::unique_ptr<std::ofstream>(
		new std::ofstream(Name));
		_jsonfs = std::unique_ptr<std::ofstream>(
		new std::ofstream(BckupName));
		
		if(_world.rank()!=0)
		{
			return;
		}

		for(int i = 0; i<numwalkers; i++)
			_root["driver"][i] = "@include(" + _prefix + "_" + std::to_string(i) + ".json)";
		_jsonfs = std::unique_ptr<std::ofstream>(
			new std::ofstream(_prefix + ".json")
		);

		_jsonfs->precision(20);

		*_jsonfs;
		Json::StreamWriterBuilder builder;
		builder["commentStyle"] = "None";
		builder["indentation"] = "\t";
		std::unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());
		writer->write(_root, &*_jsonfs);
		*_jsonfs << std::endl;

		_jsonfs->close();
		_root.clear();
	}

	void JSONObserver::PreVisit()
	{
		std::string Name = _prefix + "_" + std::to_string(_wid) +".json";
		std::string BckupName = Name + "_backup";
		std::string tempName = BckupName + "_temp";

		// Move backup to temp
		int result = rename(BckupName.c_str(), tempName.c_str());
		if(result != 0)
			throw std::runtime_error("Could not move backup file to temp file!");
		
		// Move old file to backup
		result = rename(Name.c_str(),BckupName.c_str());
		if(result != 0)
			throw std::runtime_error("Could not move restart file to backup file!");
		
		// Create pointer to New file
		_jsonfs = std::unique_ptr<std::ofstream>(
			new std::ofstream(Name));

		// Delete temp file
		result = remove(tempName.c_str());
		if(result != 0)
			throw std::runtime_error("Could not remove temp file!");

		_jsonfs->precision(20);
	}

	void JSONObserver::PostVisit()
	{
		*_jsonfs;
		Json::StreamWriterBuilder builder;
		builder["commentStyle"] = "None";
		builder["indentation"] = "\t";
		std::unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());
		writer->write(_root, &*_jsonfs);
		*_jsonfs << std::endl;

		_jsonfs->close();
		_root.clear();
	}

	void JSONObserver::Visit(const Driver& d)
	{
		d.Serialize(_root);
	}

	void JSONObserver::Serialize(Json::Value& json) const
	{
		json["type"] = "JSON";
		json["frequency"] = GetFrequency();
		json["file name"] = _prefix;
	}

}