/**
 * This file has been obtained from
 * SAPHRON - Statistical Applied PHysics through Random On-the-fly Numerics
 * https://github.com/hsidky/SAPHRON
 *
 * Copyright 2016 Hythem Sidky
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE. 
*/
#pragma once 

#include "Requirement.h"
#include "StringRequirement.h"
#include "IntegerRequirement.h"
#include "NumberRequirement.h"
#include "ObjectRequirement.h"
#include "RequirementLoader.h"
#include <algorithm>
#include <set>

namespace Json
{
	//! Array of Requirements.
	/*!
	 * \ingroup Json
	 */
	class ArrayRequirement : public Requirement
	{
	private:
		RequireList items_; //!< List of Requirements.
		std::unique_ptr<Requirement> item_; //!< Single Requirement.

		bool addItems_; //!< \c True if items can be added.
		bool setMin_; //!< \c True if minimum has been set.
		bool setMax_; //!< \c True if maximum has been set.
		bool unique_; //!< \c True if all Requirements are unique.

		int min_; //!< Minimum value.
		int max_; //!< Maximum value.

		//! Test if set is unique
		/*!
		 * \tparam T Container class.
		 * \param X Container containing Requirements.
		 * \return True if no element in the container is duplicate.
		 */
		template <class T>
		bool IsUnique(T X) 
		{
		  std::set<T> Y(X.begin(), X.end());
		  return X.size() == Y.size();
		}

	public:
		//! Constructor
		ArrayRequirement() : 
		items_(0), item_(nullptr), addItems_(true),
		setMin_(false), setMax_(false), unique_(false), 
		min_(0), max_(0) {}
		~ArrayRequirement() 
		{
			items_.clear();
		}

		//! Clear list of error messages for all Requirements.
		virtual void ClearErrors() override
		{
			if(item_ != nullptr)
				item_->ClearErrors();

			for(auto& c : items_)
				c->ClearErrors();

			Requirement::ClearErrors();
		}

		//! Clear notices for all Requirements.
		virtual void ClearNotices() override
		{
			if(item_ != nullptr)
				item_->ClearNotices();

			for(auto& c : items_)
				c->ClearNotices();

			Requirement::ClearNotices();
		}

		//! Reset array.
		virtual void Reset() override
		{
			items_.clear();

			item_.reset();

			addItems_ = true;
			setMin_ = setMax_ = unique_ = false;
			min_ = max_ = 0;
		}

		//! Parse JSON value and fill array.
		/*!
		 * \param json JSON input value.
		 * \param path Path for JSON path specification.
		 *
		 * This function parses the JSON input value and creates new
		 * Requirements for all values in the "items" branch.
		 */
		virtual void Parse(Value json, const std::string& path) override
		{
			Reset();
			RequirementLoader loader;
			
			if(json.isMember("items") && json["items"].isObject())
			{
				if((item_ = loader.LoadRequirement(json["items"])))
					item_->Parse(json["items"], path);				
			}
			else if(json.isMember("items") && json["items"].isArray())
			{
				for(auto& item : json["items"])
				{			
					if(auto iptr = loader.LoadRequirement(item))
					{
						items_.push_back(std::move(iptr));
						items_.back()->Parse(item, path);
					}
				}
			}

			if(json.isMember("additionalItems") && json["additionalItems"].isBool())
				addItems_ = json["additionalItems"].asBool();

			if(json.isMember("minItems") && json["minItems"].isInt())
			{
				setMin_ = true;
				min_ = json["minItems"].asInt();
			}

			if(json.isMember("maxItems") && json["maxItems"].isInt())
			{
				setMax_ = true;
				max_ = json["maxItems"].asInt();
			}

			if(json.isMember("uniqueItems") && json["uniqueItems"].isBool())
				unique_ = json["uniqueItems"].asBool();
		}

		//! Validate json value.
		/*!
		 * \param json JSON value to be validated.
		 * \param path Path for JSON path specification.
		 */
		virtual void Validate(const Value& json, const std::string& path) override
		{
			if(!json.isArray())
			{
				PushError(path + ": Must be of type \"array\"");
				return;
			}

			if(item_ != nullptr)
			{
				int i = 0;
				for(const auto& item : json)
				{
					item_->Validate(item, path + "/" + std::to_string(i));
					++i;
				}

				// Merge errors.
				for(const auto& error : item_->GetErrors())
					PushError(error);

				for(const auto& notice : item_->GetNotices())
					PushNotice(notice);
			}

			if(items_.size() != 0)
			{
				
				if(!addItems_ && json.size() >  items_.size())
					PushError(path + ": There cannot be any additional items in the array");

				int iv = std::min((int)items_.size(), (int)json.size());
				for(int i = 0; i < iv; ++i)
				{
					items_[i]->Validate(json[i], path);

					// Merge errors.
					for(const auto& error : items_[i]->GetErrors())
						PushError(error);

					for(const auto& notice : items_[i]->GetNotices())
						PushNotice(notice);
				}					
			}

			if(setMin_ && (int)json.size() < min_)
				PushError(path + ": There must be at least " + std::to_string(min_) + " elements in the array");
				
			if(setMax_ && (int)json.size() > max_)
				PushError(path + ": There must be no more than " + std::to_string(max_) + " elements in the array");

			if(unique_ && !IsUnique(json))
				PushError(path + ": Entries must be unique");


		}
	};
}