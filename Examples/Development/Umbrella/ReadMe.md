Method and CV Development
============
This is an example of how to develop methods and CVs in SSAGES with JSON schema.

NOTE: ALL FILES IN THIS DIRECTORY ARE FOR ILLUSTRATIVE EXAMPLE ONLY. CORRECT IMPLEMENTATION AND MODIFICATION WILL NEED TO BE DONE IN THE CORRESPONDING "schema" OR "src" DIRECTORIES.

We will be taking a look at the Umbrella sampling method and the Torsional CV.

##Universal Development
When adding a new method or CV there are a few universal steps you must take:

1. Code your method or CV
2. Code in the construction of your method or CV into either Method.cpp or CollectiveVariable.cpp.
3. Create a schema for your method or CV

##CV Development
Lets take the torsional cv as an example:

1. See TorsionalCV.h on how it is coded
2. In CollectiveVariable.cpp you need to add the JSON interpreter in the function:

```C++
CollectiveVariable* CollectiveVariable::BuildCV(const Value &json, 
					const std::string& path)
```

Using our torsional example:

```C++
else if(type == "Torsional")
{
	reader.parse(JsonSchema::TorsionalCV, schema);
	validator.Parse(schema, path);

	// Validate inputs.
	validator.Validate(json, path);
	if(validator.HasErrors())
		throw BuildException(validator.GetErrors());

	std::vector<int> atomids;
	for(auto& s : json["atom ids"])
		atomids.push_back(s.asInt());

	auto periodic = json.get("periodic", false).asBool();

	auto* c = new TorsionalCV(atomids[0], atomids[1], atomids[2], atomids[3], periodic);

	cv = static_cast<CollectiveVariable*>(c);
}
```
3. Create the required schema (see torsional.cv.json)

##Method Development
Similar steps for method development as well:

1. See Umbrella.h on how it is coded
2. In Method.cpp you need to add the JSON interpreter in the function:

```C++
Method* Method::BuildMethod(const Value &json, 
					boost::mpi::communicator& world, 
					boost::mpi::communicator& comm,
					const std::string& path)
```

Using our umbrella method example:

```C++
if(type == "Umbrella")
{
	reader.parse(JsonSchema::UmbrellaMethod, schema);
	validator.Parse(schema, path);

	// Validate inputs.
	validator.Validate(json, path);
	if(validator.HasErrors())
		throw BuildException(validator.GetErrors());

	std::vector<double> ksprings;
	for(auto& s : json["ksprings"])
		ksprings.push_back(s.asDouble());

	std::vector<double> centers;
	for(auto& s : json["centers"])
		centers.push_back(s.asDouble());

	if(ksprings.size() != centers.size())
		throw BuildException({"Need to define a spring for every center or a center for every spring!"});

	auto freq = json.get("frequency", 1).asInt();

	auto* m = new Umbrella(world, comm, ksprings, centers, freq);

	method = static_cast<Method*>(m);
}
```

3. Create the required schema (see umbrella.method.json)