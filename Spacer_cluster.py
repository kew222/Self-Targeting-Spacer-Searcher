#Spacer_cluster  Reorders spacers so that they are clustered together (if they have the same PAM

import getopt
import sys

help_message = '''
Spacer_collapse.py <file.txt>

Options:
    -h, --help    Raise this help message

*Note: Only keeps the first instance of         
'''
class Usage(Exception):
    def __init__(self,msg):
        self.msg = msg

class Params:
    def __init__(self):
        pass

    def parse_options(self, argv):
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help"])
        except getopt.error as msg:
            raise Usage(msg)
        
        for option, value in opts:
            if option in ("-h", "--help"):
                raise Usage(help_message)  
                                                                    
        if len(args) != 1:
            raise Usage(help_message)    
                                    
        return args

def main(argv=None):
    
    params = Params()     
    try:
        if argv is None:
            argv = sys.argv
            args = params.parse_options(argv)
            
            spacer_file = args[0]
            
            with open(spacer_file, 'r') as file1:
                lines = file1.readlines()
            #Change when output changes with new code!
            #Cluster the spacers based on having the same PAMs and targets
            collapsed = []; collapsed_PAMs_targets = []; num_spacers = 0; PAMs_lookup = {}
            for line in lines:
                PAMs_target = line.split('\t')[6:8]
                if PAMs_target not in collapsed_PAMs_targets:
                    PAMs_lookup["".join(PAMs_target)] = num_spacers
                    collapsed.append([line])    
                    collapsed_PAMs_targets.append(PAMs_target)
                    num_spacers += 1
                
                else:
                    collapsed[PAMs_lookup["".join(PAMs_target)]].append(line) 
                    
            with open('Clustered_spacers.txt','w') as file1:
                for x in collapsed:
                    for y in x:
                        file1.write(y)
        
    except Usage as err:
        print(sys.argv[0].split("/")[-1] + ": " + str(err.msg))
        return 2

if __name__ == "__main__":
    sys.exit(main())