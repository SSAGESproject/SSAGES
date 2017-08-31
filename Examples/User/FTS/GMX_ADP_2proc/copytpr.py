from shutil import copyfile

src = "./adp.tpr"

num_nodes = 22

for i in range(0, num_nodes):
    dest = "./adp" + str(i) + ".tpr"
    copyfile(src, dest)


