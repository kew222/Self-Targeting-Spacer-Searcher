filename = raw_input("What should the filename be? ")

with open("temp.txt", 'rU') as file1:
    lines = file1.readlines()

with open("{0}.txt".format(filename), "w") as file1:
    for line in lines:
        file1.write(line.replace('\t','\n'))
