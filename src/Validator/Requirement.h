#pragma once 

#include <vector>
#include <map>
#include <memory>
#include "json/json.h"

namespace Json
{
	class Requirement
	{
	private:
		std::vector<std::string> _errors;
		std::vector<std::string> _notices;

	protected:
		void PushError(const std::string& error) { _errors.push_back(error); }

		void PushNotice(const std::string& notice) { _notices.push_back(notice); }
	
	public:
		virtual void Parse(Value json, const std::string& path) = 0;

		virtual void Validate(const Value& json, const std::string& path) = 0;

		virtual void Reset() = 0;

		bool HasErrors() { return _errors.size() != 0; };

		std::vector<std::string> GetErrors() { return _errors; };

		virtual void ClearErrors() { _errors.clear(); }

		virtual bool HasNotices() {return _notices.size() != 0; };

		std::vector<std::string> GetNotices() {return _notices; };

		virtual void ClearNotices() { _notices.clear(); }

		virtual ~Requirement() {}
	};

	using RequireList = std::vector<std::unique_ptr<Requirement>>;
}