import matplotlib.pyplot as plt 

true_params = dict()
fitted_params = dict()
fitted_param_errors = dict()

with open("mp2test1_true_params.txt", 'r') as true_f:
    for line in true_f.readlines():
        items = line.split()
        
        try: 
            int(items[0])
            true_params[int(items[0])] = float(items[1])
        except:
            print "not a param line"



with open("millepede.res", 'r') as fitted_f:
    for line in fitted_f.readlines():
        items = line.split()
        
        try: 
            int(items[0])
            fitted_params[int(items[0])] = float(items[1])
            fitted_param_errors[int(items[0])] = float(items[4])
        except:
            print "not a param line"


errors = []

for key in sorted(true_params.iterkeys()):
    
    print key, (fitted_params[key] - true_params[key]) / fitted_param_errors[key]
