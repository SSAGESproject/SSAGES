#pragma once 

#include <iostream>

namespace SSAGES
{
	class JsonSchema
	{
	public:
		//INSERT_DEC_HERE
		static std::string grid;
		static std::string Simulation;
		static std::string driver;
		static std::string LAMMPSDriver;
		static std::string TorsionalCV;
		static std::string ImproperCV;
		static std::string AtomSeparationCV;
		static std::string AtomPositionCV;
		static std::string AtomCoordinateCV;
		static std::string CVs;
		static std::string UmbrellaMethod;
		static std::string SwarmMethod;
		static std::string methods;
		static std::string MetadynamicsMethod;
		static std::string FTSMethod;
		static std::string ForwardFluxMethod;
		static std::string ElasticBandMethod;
		static std::string BFSMethod;
		static std::string ABFMethod;
		
	};
}