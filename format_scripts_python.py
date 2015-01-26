


########################## 
############
######
# python format_scripts_python.py dreme2meme('dreme.txt')
######
def dreme2meme(name):
    data0 = open(name,'r')

    data1 = []

    for records in data0:
        if records.split()!=[]:
            if records.split()[0] != '#':
                data1.append(records)
        else:
            data1.append(records)
    result = open(name+'.meme','w')

    for records in data1:
        result.write(records)

    result.close()
    data0.close()
############
##########################





