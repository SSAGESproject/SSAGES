#pragma once 

#include "Requirement.h"
#include "RequirementLoader.h"

namespace Json
{
	class NotRequirement : public Requirement
	{
	private:
		std::unique_ptr<Requirement> _req;

	public:
		NotRequirement() : _req(nullptr) {}

		virtual void ClearErrors() override
		{
			if(_req)
				_req->ClearErrors();

			Requirement::ClearErrors();
		}

		virtual void ClearNotices() override
		{
			if(_req != nullptr)
				_req->ClearNotices();
		
			Requirement::ClearNotices();
		} 

		virtual void Reset() override
		{
			ClearErrors();
			ClearNotices();
			
			_req.reset();		
		}

		virtual void Parse(Value json, const std::string& path) override
		{
			Reset();
			RequirementLoader loader;

			auto& head = json.isMember("not") ? json["not"] : json;
			if((_req = loader.LoadRequirement(head)))
				_req->Parse(head, path);
		}

		virtual void Validate(const Value& json, const std::string& path) override
		{
			_req->Validate(json, path);
			if(!_req->HasErrors())
				PushError(path + ": Value must not validate against requirement.");
		}
	};
}