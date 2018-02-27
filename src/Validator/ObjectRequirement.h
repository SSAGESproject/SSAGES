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

#include <iostream>
#include <algorithm>
#include <list>
#include <map>
#include <numeric>
#include "Requirement.h"
#include "StringRequirement.h"
#include "IntegerRequirement.h"
#include "NumberRequirement.h"
#include "DependencyRequirement.h"
#include "RequirementLoader.h"

namespace Json
{
	//! Requirements on an object
	/*!
	 * \ingroup Json
	 */
	class ObjectRequirement : public Requirement
	{
	private:
		//! Map of properties the object needs to have.
		std::map<std::string, std::unique_ptr<Requirement>> properties_;

		//! Map of patterns the object needs to match.
		std::map<std::string, std::unique_ptr<Requirement>> patternProps_;

		//! List of requirements.
		RequireList extended_;

		//! Dependency requirement.
		std::unique_ptr<DependencyRequirement> dependency_;

		std::vector<std::string> required_; //!< List of requirements.
		bool moreProps_; //!< If \c True, more properties need to be set.
		bool setMin_; //!< If \c True lower bound is active.
		bool setMax_; //!< If \c True upper bound is active.
		unsigned int min_; //!< Lower bound.
		unsigned int max_; //!< Upper bound.

	public:
		//! Constructor.
		ObjectRequirement() : 
		properties_(), patternProps_(), extended_(0), dependency_(nullptr), required_(),
		moreProps_(true), setMin_(false), setMax_(false), min_(0), max_(0)
		{}

		//! Destructor.
		~ObjectRequirement()
		{
			properties_.clear();

			patternProps_.clear();

			extended_.clear();
		}

		//! Clear errors on all Requirements.
		virtual void ClearErrors() override
		{
			for(auto& c : properties_)
				c.second->ClearErrors();

			for(auto& c : patternProps_)
				c.second->ClearErrors();

			for(auto& c : extended_)
				c->ClearErrors();

			if(dependency_ != nullptr)
				dependency_->ClearErrors();

			Requirement::ClearErrors();
		}

		//! Clear notices on all Requirements.
		virtual void ClearNotices() override
		{
			for(auto& c : properties_)
				c.second->ClearNotices();

			for(auto& c : patternProps_)
				c.second->ClearNotices();

			for(auto& c : extended_)
				c->ClearNotices();

			if(dependency_ != nullptr)
				dependency_->ClearNotices();

			Requirement::ClearNotices();
		} 

		//! Reset Requirement.
		virtual void Reset() override
		{
			ClearErrors();
			ClearNotices();

			properties_.clear();
			patternProps_.clear();

			moreProps_ = true;
			setMin_ = setMax_ = false;
			min_ = max_ = 0;
			required_.clear();
			dependency_.reset();
		}

		//! Parse JSON value to generate Requirement(s).
		/*!
		 * \param json JSON input value.
		 * \param path Path for JSON path specification.
		 */
		virtual void Parse(Value json, const std::string& path) override
		{
			Reset();
			RequirementLoader loader;

			// Additional properties.
			if(json.isMember("additionalProperties"))
			{
				if(json["additionalProperties"].isBool())
					moreProps_ = json["additionalProperties"].asBool();
				else if(json["additionalProperties"].isObject())
					json["properties"]["additionalProperties"] = json["additionalProperties"];
			}

			// properties.
			if(json.isMember("properties") && json["properties"].isObject())
			{
				auto& props = json["properties"];
				auto names = props.getMemberNames();
				int i = 0;
				for(auto& prop : props)
				{
					if(auto property = loader.LoadRequirement(prop))
					{
						properties_[names[i]] = std::move(property);
						properties_[names[i]]->Parse(prop, path + "/" + names[i]);
					}

					++i;
				}
			}

			// Pattern properties. TODO: eliminate redundant code!!
			if(json.isMember("patternProperties") && json["patternProperties"].isObject())
			{
				auto& props = json["patternProperties"];
				auto names = props.getMemberNames();
				int i = 0;
				for(auto& prop : props)
				{
					if(prop.isObject())
					{
						if(auto property = loader.LoadRequirement(prop))
						{
							patternProps_[names[i]] = std::move(property);
							patternProps_[names[i]]->Parse(prop, path + "/" + names[i]);
						}
					}

					++i;
				}
			}

			// Required properties.
			if(json.isMember("required") && json["required"].isArray())
			{
				for(auto& requirement : json["required"])
					required_.push_back(requirement.asString());
			}

			// Min property count.
			if(json.isMember("minProperties") && json["minProperties"].isUInt())
			{
				setMin_ = true;
				min_ = json["minProperties"].asInt();
			}

			// Max property count.
			if(json.isMember("maxProperties") && json["maxProperties"].isUInt())
			{
				setMax_ = true;
				max_ = json["maxProperties"].asInt();
			}

			// Dependencies
			if(json.isMember("dependencies") && json["dependencies"].isObject())
			{
				dependency_ = std::unique_ptr<DependencyRequirement>(new DependencyRequirement());
				dependency_->Parse(json["dependencies"], path);
			}

			// Extended properties. 
			for(auto& prop : json)
			{
				if(auto req  = loader.LoadExtended(prop))
				{
					extended_.push_back(std::move(req));
					extended_.back()->Parse(prop, path);
				}
			}
		}

		//! Validate JSON value.
		/*!
		 * \param json JSON value to be validated.
		 * \param path Path for JSON path specification.
		 */
		virtual void Validate(const Value& json, const std::string& path) override
		{
			if(!json.isObject())
			{
				PushError(path + ": Must be of type \"object\"");
				return;
			}

			if(setMin_ && json.size() < min_)
				PushError(path + ": Object must contain at least " + std::to_string(min_) + " properties");

			if(setMax_ && json.size() > max_)
				PushError(path + ": Object must contain at most " + std::to_string(max_) + " properties");


			// Check dependencies. 
			if(dependency_ != nullptr)
			{
				dependency_->Validate(json, path);
				if(dependency_->HasErrors())
					for(const auto& error : dependency_->GetErrors())
						PushError(error);
				if(dependency_->HasNotices())
					for(const auto& notice : dependency_->GetNotices())
						PushNotice(notice);
			}

			// Copy so we can pop items off the list.
			auto rprops = required_;

			auto names = json.getMemberNames();
			int i = 0;
			for(auto& prop : json)
			{
				Requirement* requirement = nullptr; 
				auto it = properties_.find(names[i]);
				if(it != properties_.end())
					requirement = it->second.get();					
				else if(patternProps_.size() != 0)
				{
					for(auto& pattern : patternProps_)
					{
						auto regex = std::regex(pattern.first, std::regex::ECMAScript);
						if(std::regex_search(names[i], regex))
							requirement = pattern.second.get();
					}
				}
				
				if(!requirement && properties_.find("additionalProperties") != properties_.end())
					requirement = properties_["additionalProperties"].get();
				
				if(requirement)
				{
					requirement->Validate(prop, path + "/" + names[i]);
					if(requirement->HasErrors())
						for(const auto& error : requirement->GetErrors())
							PushError(error);
					if(requirement->HasNotices())
						for(const auto& notice : requirement->GetNotices())
							PushNotice(notice);
				}
				else if(!moreProps_)
					PushError(path + ": Invalid property \"" + names[i] + "\" specified");

				rprops.erase(std::remove(rprops.begin(), rprops.end(),names[i]),rprops.end());
				++i;
			}

			if(required_.size() && rprops.size() != 0)
			{
				std::string msg = std::accumulate(rprops.begin(), rprops.end(), std::string(), 
    				[](const std::string& a, const std::string& b) -> std::string { 
        				return a + (a.length() > 0 ? ", " : "") + b; 
    			});
				PushError(path + ": Missing properties: " + msg);
			}

			// Validate extended.
			for(auto& requirement : extended_)
			{
				requirement->Validate(json, path);
				if(requirement->HasErrors())
						for(const auto& error : requirement->GetErrors())
							PushError(error);
				if(requirement->HasNotices())
					for(const auto& notice : requirement->GetNotices())
						PushNotice(notice);
			}
		}
	};
}
