def import_data():
	input_file = "Spacers_outside_islands.txt"
	with open(input_file, 'r') as fileread:
			lines = fileread.readlines()
	imported_data = []
	for line in lines:
		imported_data.append(line.strip().split('\t'))


def reverse_complement(sequence):
	sequence = imported_data(0)(7)
	reverse_complement= sequence[::-1]
	str = imported_data(0)(7)
	print str.replace("T","A")
	print str.replace("A","T")
	print str.replace("G","C")
	print str.replace("C","G")
	reverse_complement= sequence[::-1]
	print reverse_complement
	
	


