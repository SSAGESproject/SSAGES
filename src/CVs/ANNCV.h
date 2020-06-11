/**
 * This file is part of
 * SSAGES - Software Suite for Advanced General Ensemble Simulations
 *
 * Copyright 2020 Wei Chen <wei.herbert.chen@gmail.com>
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

#include "CollectiveVariable.h"
#include "Validator/ObjectRequirement.h"
#include "Drivers/DriverException.h"
#include "Snapshot.h"
#include "schema.h"

#include <array>
#include <cmath>
#include <assert.h>

namespace SSAGES
{
	//! ANN (artifical neural network) collective variables.
	/*! This CV takes scaled (specified by "scaling_factor")
	* Cartesian coordinates of a group of atoms (specified by "atomids") as inputs to a neural network
	* (its number of nodes, connection weights, and activation functions are specified by "num_nodes",
	* "coeff_file", "activations", respectively),
	* computes one component (specified by "out_index") of the final neural network outputs as the CV value.
	*
	*/
	class ANNCV : public CollectiveVariable
	{
	private:
		Label atomids_;   // indices of atoms
		double scaling_factor_;  // scaling factor to make input unitless, note that it is unit dependent (a model trained with A will have different scaling factor with one trained with nm)
		std::vector<int> num_nodes_;  // numbers of nodes for neural network
		std::vector<Eigen::MatrixXd> weight_coeff_;
		std::vector<Vector> bias_;
		std::vector<std::string> activations_;
		int out_index_;     // index of output component

	public:

		//! Constructor.
		/*!
		 * \param atomids atom IDs used as inputs to ANN
		 * \param scaling_factor scaling factor for input coordinates (input coordinates will be divided by this factor)
		 * \param num_nodes  numbers of nodes for each layer of ANN
		 * \param coeff_file file storing weights and bias of ANN
		 * \param activations activations functions (as strings) of ANN
		 * \param out_index index of output component of ANN
		 *
		 * \todo Bounds needs to be an input and periodic boundary conditions
		 */
		ANNCV(
			Label atomids,
			double scaling_factor,
			std::vector<int> num_nodes,
			std::string coeff_file, // file storing weights and bias of ANN
			std::vector<std::string> activations,
			int out_index
		) :
			atomids_(atomids), scaling_factor_(scaling_factor), num_nodes_(num_nodes),
			activations_(activations), out_index_(out_index)
	  	{
			// read coefficients from file
			std::ifstream my_f(coeff_file);
			std::string temp_vec;
			int layer_index = 0;
			if (num_nodes_[0] != atomids_.size() * 3) {
				throw BuildException({
					"WARNING: input dim should be " + std::to_string(atomids_.size() * 3) + " found: "
					+ std::to_string(num_nodes_[0])
				});
			}
			while (std::getline(my_f, temp_vec) )
			{
				std::istringstream ss(temp_vec);
				std::string token;
				std::vector<double> temp_weight, temp_bias;
				while (std::getline(ss, token, ','))  // coefficients are separated by comma
				{
					temp_weight.push_back(stod(token));
				}
				if (temp_weight.size() != num_nodes_[layer_index] * num_nodes_[layer_index + 1]) {
					throw BuildException({
						"WARNING: layer weight size = " + std::to_string(temp_weight.size()) +  " expected: "
						+ std::to_string(num_nodes_[layer_index] * num_nodes_[layer_index + 1])
					});
				}
				Eigen::Map<Eigen::MatrixXd> temp_weight_v(
					&temp_weight[0], num_nodes_[layer_index], num_nodes_[layer_index + 1]
				);
				std::getline(my_f, temp_vec);
				std::istringstream ss2(temp_vec);
				while (std::getline(ss2, token, ','))
				{
					temp_bias.push_back(stod(token));
				}
				if (temp_bias.size() != num_nodes_[layer_index + 1]) {
					throw BuildException({
						"WARNING: layer bias size = " + std::to_string(temp_bias.size())
						+ " expected: " + std::to_string(num_nodes_[layer_index + 1])
					});
				}
				Vector temp_bias_v = Vector::Map(temp_bias.data(), temp_bias.size());
				weight_coeff_.push_back(temp_weight_v);
				bias_.push_back(temp_bias_v);
				layer_index ++;
			}
	  	}

		//! Initialize necessary variables.
		/*!
		 * \param snapshot Current simulation snapshot.
		 */
		void Initialize(const Snapshot& snapshot) override
		{
			using std::to_string;

			std::vector<int> found;
			snapshot.GetLocalIndices(atomids_, &found);
			int nfound = found.size();
			MPI_Allreduce(MPI_IN_PLACE, &nfound, 1, MPI_INT, MPI_SUM, snapshot.GetCommunicator());

			if(nfound != atomids_.size())
				throw BuildException({
					"ANNCV: Expected to find " +
					to_string(atomids_.size()) +
					" atoms, but only found " +
					to_string(nfound) + "."
				});
		}

		std::vector<Vector> forward_prop(Vector& input_vec) {
			std::vector<Vector> output_of_layers;
			Vector temp_out = input_vec;
			output_of_layers.push_back(temp_out);
			for (int ii = 0; ii < weight_coeff_.size(); ii ++) {
				temp_out = weight_coeff_[ii].transpose() * temp_out + bias_[ii];
				if (activations_[ii] == "Tanh") {
					for (int kk = 0; kk < temp_out.size(); kk ++) {
						temp_out[kk] = std::tanh(temp_out[kk]);
					}
				}
				else if (activations_[ii] == "ReLU") {
					for (int kk = 0; kk < temp_out.size(); kk ++) {
						temp_out[kk] = temp_out[kk] < 0 ? 0 : temp_out[kk];
					}
				}
				output_of_layers.push_back(temp_out);
			}
			return output_of_layers;
		}

		std::vector<Vector> back_prop(std::vector<Vector>& output_of_layers) {
			auto deriv_back = output_of_layers;
			int num = output_of_layers.size();
			for (int ii = 0; ii < output_of_layers[num - 1].size(); ii ++ ) {
				if (ii == out_index_) {
					deriv_back[num - 1][ii] = 1;
				}
				else {
					deriv_back[num - 1][ii] = 0;
				}
			}
			for (int ii = num - 2; ii >= 0; ii --) {
				if (activations_[ii] == "Tanh") {
					for (int kk = 0; kk < deriv_back[ii + 1].size(); kk ++) {
						deriv_back[ii + 1][kk] = deriv_back[ii + 1][kk] * (
							1 - output_of_layers[ii + 1][kk] * output_of_layers[ii + 1][kk]);
					}
				}
				deriv_back[ii] = weight_coeff_[ii] * deriv_back[ii + 1];
			}
			return deriv_back;
		}

		//! Evaluate the CV.
		/*!
		 * \param snapshot Current simulation snapshot.
		 */
		void Evaluate(const Snapshot& snapshot) override
		{
			// Get data from snapshot.
			auto n = snapshot.GetNumAtoms();
			const auto& pos = snapshot.GetPositions();
			auto& comm = snapshot.GetCommunicator();

			// Initialize gradient.
			std::fill(grad_.begin(), grad_.end(), Vector3{0,0,0});
			grad_.resize(n, Vector3{0,0,0});

			// Vector3 xi{0, 0, 0}, xj{0, 0, 0}, xk{0, 0, 0};
			std::vector<Vector3> positions(atomids_.size(), Vector3({0,0,0}));
			Label local_idx;
			snapshot.GetLocalIndices(atomids_, &local_idx);
			for (int ii = 0; ii < atomids_.size(); ii ++) {
				if (local_idx[ii] != -1) {
					positions[ii] = pos[local_idx[ii]];
				}
			}
			// By performing a reduce, we actually collect all. This can
			// be converted to a more intelligent allgater on rank then bcast.
			MPI_Allreduce(MPI_IN_PLACE, positions.data(), positions.size(), MPI_DOUBLE, MPI_SUM, comm);
			auto com = snapshot.CenterOfMass(atomids_, false);  // center of mass coordinates (not weighted with mass)
			// remove translation degree of freedom
			for (auto& temp_pos: positions) {
				temp_pos = temp_pos - com;
			}
			Vector input_vec(positions.size() * 3);  // flatten input vector
			for (int ii = 0; ii < positions.size(); ii ++) {
				for (int jj = 0; jj < 3; jj ++) {
					input_vec[ii * 3 + jj] = positions[ii][jj];
				}
			}
			input_vec = input_vec / scaling_factor_;
			auto output_of_layers = forward_prop(input_vec);
			val_ = output_of_layers[output_of_layers.size() - 1][out_index_];
			auto deriv_back = back_prop(output_of_layers);
			// subtract mean from deriv_back[0]
			int num_atoms = deriv_back[0].size() / 3;
			double average[3] = {0};
			for (int kk = 0; kk < 3; kk ++) {
				for (int ii = 0; ii < num_atoms; ii ++) {
					average[kk] += deriv_back[0][ii * 3 + kk];
				}
				average[kk] /= num_atoms;
			}
		}

		//! \copydoc CollectiveVariable::BuildCV()
		static ANNCV* Build(const Json::Value& json, const std::string& path)
		{
			Json::ObjectRequirement validator;
			Json::Value schema;
			Json::CharReaderBuilder rbuilder;
			Json::CharReader* reader = rbuilder.newCharReader();

			reader->parse(JsonSchema::ANNCV.c_str(),
			              JsonSchema::ANNCV.c_str() + JsonSchema::ANNCV.size(),
			              &schema, NULL);
			validator.Parse(schema, path);

			// Validate inputs.
			validator.Validate(json, path);
			if(validator.HasErrors())
				throw BuildException(validator.GetErrors());

			std::vector<int> atomids;
			for(auto& s : json["atom_ids"])
				atomids.push_back(s.asInt());
			double scaling_factor = json["scaling_factor"].asDouble();
			std::vector<int> num_nodes;
			for (auto &s : json["num_nodes"])
				num_nodes.push_back(s.asInt());
			std::vector<std::string> activations;
			std::string coeff_file = json["coeff_file"].asString();
			for (auto &s : json["activations"]) {
				auto name = s.asString();
				if (name == "Tanh" || name == "ReLU" || name == "Linear")
					activations.push_back(name);
				else
					throw std::runtime_error("invalid activation function " + name + " provided");
			}
			int out_index = json["index"].asInt();
			return new ANNCV(atomids, scaling_factor, num_nodes, coeff_file, activations, out_index);
		}
	};
}
