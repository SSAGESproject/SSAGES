#include <stdexcept>
#include "Validator/ObjectRequirement.h"
#include "Drivers/DriverException.h"
#include "Basis.h"
#include "schema.h"

namespace SSAGES
{
    BasisFunction* BasisFunction::Build(const Json::Value& json, const std::string& path, const int nbins)
    {
        auto type = json.get("type","none").asString();
        if(type == "Legendre")
            return Legendre::Build(json, path, nbins);
        else if (type == "Chebyshev")
            return Chebyshev::Build(json, path, nbins);
        else
            throw std::invalid_argument("Invalid basis set type \"" + type + "\".");
    }

    Chebyshev* Chebyshev::Build(const Json::Value& json,
                              const std::string& path,
                              const int nbin)
    {
        Json::ObjectRequirement validator;
        Json::Value schema;
        Json::Reader reader;
        
        reader.parse(JsonSchema::ChebyshevBasis, schema);
        validator.Parse(schema, path);

        //Validate Inputs
        validator.Validate(json, path);
        if(validator.HasErrors())
            throw BuildException(validator.GetErrors());

        return new Chebyshev(
                json["polynomial_order"].asInt(),
                json["lower_bound"].asDouble(),
                json["upper_bound"].asDouble(),
                nbin
            );
    }

    Legendre* Legendre::Build(const Json::Value& json,
                              const std::string& path,
                              const int nbin)
    {
        Json::ObjectRequirement validator;
        Json::Value schema;
        Json::Reader reader;
        
        reader.parse(JsonSchema::ChebyshevBasis, schema);
        validator.Parse(schema, path);

        //Validate Inputs
        validator.Validate(json, path);
        if(validator.HasErrors())
            throw BuildException(validator.GetErrors());

        return new Legendre(
                json["polynomial_order"].asInt(),
                nbin
            );
    }
};
