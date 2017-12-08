import argparse

parser = argparse.ArgumentParser(description='This script clean a .xml output from Qbox.')
parser.add_argument('-i','--input',help='name of the input .xml file containing the trajectory to be cleaned. ',required=True)
parser.add_argument('-o','--output',help='name of the cleaned .xml file containing the trajectory. ',required=True)

args=parser.parse_args()

fout = open(args.output,'w')

fout.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?> \n")
fout.write("<fpmd:simulation xmlns:fpmd=\"http://www.quantum-simulation.org/ns/fpmd/fpmd-1.0\"> \n")


for line in open(args.input,'r'):
	if "<?xml" not in line and "<fpmd:" not in line and "[qbox]" not in line and "End of" not in line and "</fpmd" not in line:
		fout.write(line)

fout.write("[qbox] \n")
fout.write(" End of command stream \n")
fout.write("</fpmd:simulation>")

fout.close()
