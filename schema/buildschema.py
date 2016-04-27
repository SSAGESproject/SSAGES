#!/usr/bin/env python
'''
buildschema.py - Compiles individual schema files into a single complete
                 JSON schema and creates a C++ header file containing
                 the escaped schemas.

                 This is a modified version originally from github hsidky/SAPHRON
                 https://github.com/hsidky/SAPHRON.git
'''
import re
import os
from os import path
import json

# Definitions
currd = path.dirname(path.abspath(__file__))
schemah = '../include/schema.h'
schemac = '../src/JSON/schema.cpp'

folders = [
	"Methods/",
	"CVs/",
	"Drivers/",
	"Simulations/",
]

exclude = [

]


def processfile(filename):
	fpath = path.dirname(filename)
	with open(filename, 'r') as f:
		text = f.read()
		for match in re.findall('"@file[(](.*)[)]"', text):
			text = text.replace(
				'"@file(' + match + ')"',
				processfile(fpath + "/" + match)
			)
		return text


def headertimestamp():
	f = path.join(currd, schemah)
	if not path.isfile(f):
		return 0
	return path.getmtime(f)


def jsontimestamps(folders):
	maxts = 0
	for folder in folders:
		for f in os.listdir(path.join(currd, folder)):
			if f not in exclude:
				fpath = path.join(currd, folder, f)
				tstamp = path.getmtime(fpath)
				if tstamp > maxts:
					maxts = tstamp
	return maxts


def genfiles(folders):
	# Get timestamp of header file
	htime = headertimestamp()
	jtime = jsontimestamps(folders)

	# If no files have not been modded since header,
	# don't do anything.
	if htime >= jtime:
		return

	# Delete old header/cpp files if they exist.
	if path.isfile(path.join(currd, schemah)):
		os.remove(path.join(currd, schemah))
	if path.isfile(path.join(currd, schemac)):
		os.remove(path.join(currd, schemac))

	fheader = path.join(currd, "headertemplate.h")
	fcpp = path.join(currd, "sourcetemplate.cpp")

	hlines = []
	cpplines = []
	with open(fheader, "r") as f:
		hlines = f.readlines()
	with open(fcpp, "r") as f:
		cpplines = f.readlines()

	# Go through folders, load json and write to header.
	for folder in folders:
		for f in os.listdir(path.join(currd, folder)):
			if f not in exclude:
				fpath = path.join(currd, folder, f)
				try:
					jobj = json.loads(processfile(fpath))
					varname = jobj["varname"]
					jobj.pop("varname", None)
				except Exception as e:
					print("Error parsing {0}: {1}".format(fpath, e))

				# Add lines to header and cpp
				decl = '\t\tstatic std::string ' + varname + ';\n'
				i = next(j for j, s in enumerate(hlines) if 'INSERT_DEC_HERE' in s)
				hlines.insert(i+1, decl)

				defn = 'std::string SSAGES::JsonSchema::' + varname + ' = '
				i = next(j for j, s in enumerate(cpplines) if 'INSERT_DEF_HERE' in s)
				cpplines.insert(
					i+1,
					"\t" + defn + '"' + json.dumps(jobj).replace('"', '\\"') + '";\n'
				)

	with open(path.join(currd, schemah), "w") as f:
		f.writelines(hlines)

	with open(path.join(currd, schemac), "w") as f:
		f.writelines(cpplines)

genfiles(folders)
