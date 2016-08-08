#pragma once

#include "Method.h"
#include "../CVs/CollectiveVariable.h"
#include <fstream>
#include <vector>


namespace SSAGES
{
	//! Multidimensional hill
	/*!
	 * Structure representing a multidimensional hill (Gaussian) which is
	 * centered at "center" with widths "width" of height "height". A
	 * multidimensional Gaussian has one height but n centers and widths.
	 */
	struct Hill 
	{
		//! Hill center.
		std::vector<double> center;

		//! Hill width.
		std::vector<double> width;

		//! Hill height.
		double height;
		
		//! Constructs a multidimensional Hill (Gaussian)
		/*!
		 * \param center Hill center.
		 * \param sigma Hill width.
		 * \param height Hill height.
		 */
		Hill(const std::vector<double>& center, 
			 const std::vector<double>& sigma, 
			 double height) :
		 center(center), width(sigma), height(height)
		{}
	};

	//! "Vanilla" multi-dimensional Metadynamics.
	/*!
	 * Implementation of a "vanilla" multi-dimensional Metadynamics method with
	 * no bells and whistles.
	 *
	 * \ingroup Methods
	 */
	class Meta : public Method
	{
	private:	
		//! Hills.
		std::vector<Hill> _hills;

		//! Hill height.
		double _height;

		//! Hill widths.
		std::vector<double> _widths;

		//! Gridding flag and grid
		bool _isgrid;

		//! Derivatives.
		std::vector<double> _derivatives;

		//! Bias magnitude.
		double _bias;

		//! Frequency of new hills
		unsigned int _hillfreq;

		//! Adds a new hill.
		/*!
		 * \param cvs List of CVs.
		 */
		void AddHill(const CVList& cvs);

		//! Computes the bias force.
		/*!
		 * \param cvs List of CVs.
		 */
		void CalcBiasForce(const CVList& cvs);

		//! Prints the new hill to file
		/*!
		 * \param hill Hill to be printed.
		 */
		void PrintHill(const Hill& hill);

		//! Output stream for hill data.
		std::ofstream _hillsout;

	public:
		//! Constructor
		/*!
		 * \param world MPI global communicator.
		 * \param comm MPI local communicator.
		 * \param height Height of the hills to be deposited.
		 * \param widths Width of the hills to be deposited.
		 * \param hillfreq Frequency of depositing hills.
		 * \param frequency Frequency of invoking this method.
		 *
		 * Constructs an instance of Metadynamics method.
		 * \note The size of "widths" should be commensurate with the number of
		 *       CV's expected.
		 */
		Meta(boost::mpi::communicator& world,
			 boost::mpi::communicator& comm,
			 double height, 
			 const std::vector<double>& widths, 
			 unsigned int hillfreq, 
        		 unsigned int frequency,
		         bool isgrid ) : 
		Method(frequency, world, comm), _hills(), _height(height), _widths(widths), 
		  _derivatives(0), _bias(0), _hillfreq(hillfreq), _isgrid(isgrid)
		{
		}

		//! Pre-simulation hook.
		/*!
		 * \param snapshot Current simulation snapshot.
		 * \param cvs List of CVs.
		 */
		void PreSimulation(Snapshot* snapshot, const CVList& cvs) override;

		//! Post-integration hook.
		/*!
		 * \param snapshot Current simulation snapshot.
		 * \param cvs List of CVs.
		 */
		void PostIntegration(Snapshot* snapshot, const CVList& cvs) override;

		//! Post-simulation hook.
		/*!
		 * \param snapshot Current simulation snapshot.
		 * \param cvs List of CVs.
		 */
		void PostSimulation(Snapshot* snapshot, const CVList& cvs) override;

		//! \copydoc Serializable::Serialize()
		/*!
		 * \warning Serialization not implemented yet!
		 */
		void Serialize(Json::Value& json) const override
		{

		}

		//! Destructor.
		~Meta() {}
	};
}
			
