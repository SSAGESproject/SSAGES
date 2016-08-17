#pragma once 

#include "CollectiveVariable.h"

#include <array>
#include <cmath>

namespace SSAGES
{
	// Collective variable to calculate radius of gyration.

	class MockCV : public CollectiveVariable
	{
	private:

		// Current value of the CV.
		double _val;

		// Gradient of the CV, dRg/dxi.
		std::vector<Vector3> _grad;

		Vector3 _user_defined_grad;

		// Bounds on CV.
		std::array<double, 2> _bounds;

	public:

		//! Constructor.
		/*!
		 * \param periodic If \c True consider periodic boundary conditions.
		 *
		 * Construct a mock CV
		 */
		MockCV(double value, Vector3 grad, double upper, double lower) :
		_val(value), _grad(), _user_defined_grad(grad), _bounds{{upper, lower}}
		{

		}

		//! Initialize necessary variables.
		/*!
		 * \param snapshot Current simulation snapshot.
		 */
		void Initialize(const Snapshot& snapshot) override
		{
			// Initialize gradient
			auto n = snapshot.GetPositions().size();
			_grad.resize(n, _user_defined_grad);

		}

		//! Evaluate the CV.
		/*!
		 * \param snapshot Current simulation snapshot.
		 */
		void Evaluate(const Snapshot& snapshot) override
		{

		}

		//! Return the value of the CV.
		/*!
		 * \return Current value of the CV.
		 */
		double GetValue() const override 
		{ 
			return _val; 
		}

		//! Return value taking periodic boundary conditions into account.
		/*!
		 * \param Location Input value.
		 * \return Original value.
		 *
		 * The atom position CV does not take periodic boundary conditions into
		 * account. Thus, this function will always return the unmodified input
		 * value.
		 */
		double GetPeriodicValue(double Location) const override
		{
			return Location;
		}

		//! Return the gradient of the CV.
		/*!
		 * \return Gradient of the CV.
		 */
		const std::vector<Vector3>& GetGradient() const override
		{
			return _grad;
		}

		//! Return the boundaries of the CV.
		/*!
		 * \return Values of the lower and upper boundaries of this CV.
		 */
		const std::array<double, 2>& GetBoundaries() const override
		{
			return _bounds;
		}

		//! Get difference taking periodic boundary conditions into account.
		/*!
		 * \param Location Input value.
		 * \return Simple difference between current value of CV and Location.
		 *
		 * The atom position CV does not take periodic boundary conditions into
		 * account. Thus, this function returns the simple difference between
		 * the current value of the CV and the input value.
		 */
		double GetDifference(const double Location) const override
		{
			return _val - Location;
		}

		//! Serialize this CV for restart purposes.
		/*!
		 * \param json JSON value
		 */
		virtual void Serialize(Json::Value& json) const override
		{

		}

	};
}
