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
		RequireList _items; //!< List of Requirements.
		std::unique_ptr<Requirement> _item; //!< Single Requirement.

		bool _addItems; //!< \c True if items can be added.
		bool _setMin; //!< \c True if minimum has been set.
		bool _setMax; //!< \c True if maximum has been set.
		bool _unique; //!< \c True if all Requirements are unique.

		int _min; //!< Minimum value.
		int _max; //!< Maximum value.

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
		_items(0), _item(nullptr), _addItems(true),
		_setMin(false), _setMax(false), _unique(false), 
		_min(0), _max(0) {}
		~ArrayRequirement() 
		{
			_items.clear();
		}

		//! Clear list of error messages for all Requirements.
		virtual void ClearErrors() override
		{
			if(_item != nullptr)
				_item->ClearErrors();

			for(auto& c : _items)
				c->ClearErrors();

			Requirement::ClearErrors();
		}

		//! Clear notices for all Requirements.
		virtual void ClearNotices() override
		{
			if(_item != nullptr)
				_item->ClearNotices();

			for(auto& c : _items)
				c->ClearNotices();

			Requirement::ClearNotices();
		}

		//! Reset array.
		virtual void Reset() override
		{
			_items.clear();

			_item.reset();

			_addItems = true;
			_setMin = _setMax = _unique = false;
			_min = _max = 0;
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
				if((_item = loader.LoadRequirement(json["items"])))
					_item->Parse(json["items"], path);				
			}
			else if(json.isMember("items") && json["items"].isArray())
			{
				for(auto& item : json["items"])
				{			
					if(auto iptr = loader.LoadRequirement(item))
					{
						_items.push_back(std::move(iptr));
						_items.back()->Parse(item, path);
					}
				}
			}

			if(json.isMember("additionalItems") && json["additionalItems"].isBool())
				_addItems = json["additionalItems"].asBool();

			if(json.isMember("minItems") && json["minItems"].isInt())
			{
				_setMin = true;
				_min = json["minItems"].asInt();
			}

			if(json.isMember("maxItems") && json["maxItems"].isInt())
			{
				_setMax = true;
				_max = json["maxItems"].asInt();
			}

			if(json.isMember("uniqueItems") && json["uniqueItems"].isBool())
				_unique = json["uniqueItems"].asBool();
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

			if(_item != nullptr)
			{
				int i = 0;
				for(const auto& item : json)
				{
					_item->Validate(item, path + "/" + std::to_string(i));
					++i;
				}

				// Merge errors.
				for(const auto& error : _item->GetErrors())
					PushError(error);

				for(const auto& notice : _item->GetNotices())
					PushNotice(notice);
			}

			if(_items.size() != 0)
			{
				
				if(!_addItems && json.size() >  _items.size())
					PushError(path + ": There cannot be any additional items in the array");

				int iv = std::min((int)_items.size(), (int)json.size());
				for(int i = 0; i < iv; ++i)
				{
					_items[i]->Validate(json[i], path);

					// Merge errors.
					for(const auto& error : _items[i]->GetErrors())
						PushError(error);

					for(const auto& notice : _items[i]->GetNotices())
						PushNotice(notice);
				}					
			}

			if(_setMin && (int)json.size() < _min)
				PushError(path + ": There must be at least " + std::to_string(_min) + " elements in the array");
				
			if(_setMax && (int)json.size() > _max)
				PushError(path + ": There must be no more than " + std::to_string(_max) + " elements in the array");

			if(_unique && !IsUnique(json))
				PushError(path + ": Entries must be unique");


		}
	};
}