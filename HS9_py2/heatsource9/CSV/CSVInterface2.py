"""This is an alternative CSV interface that uses 'csv' instead of 'pandas'. it's not working right now."""

from __future__ import division, print_function
from os.path import exists, join, normpath
from collections import defaultdict
import csv
import pandas as pd

def csv_reader(inputdir,filename):
    with open(join(inputdir,filename), "rU") as file_object:
        data = [row for row in csv.reader(file_object.read().splitlines(),dialect="excel")]
    return data

def csv_dictreader(inputdir,filename,colnames):
    data = defaultdict(list) # each value in each column is appended to a list
    with open(join(inputdir,filename), "rU") as file_object:
        reader = csv.DictReader(file_object,dialect="excel")
        reader.next() # skip the header row
        for row in reader: # read a row as {column1: value1, column2: value2,...}
                for (k,v) in row.items(): # go over each column name and value 
                    data[k].append(v) # append the value into the appropriate list based on column name k      
    return data


############################
# Landcover Codes #
inputdir = "/Users/rmichie/Desktop/2003_PRE_CCC/heatsource_9/inputfiles/"
filename = "LC_codes_LAI.csv"
colnames = ['name', 'code','height','density','k','overhang']

#data = []
#data = csv_dictreader(inputdir, filename, colnames)

#codes = list(data['code'])
# make a list of lists with values: [(height[0], dens[0], k[0], over[0]), (height[1],...),...]
#vals = [tuple([j for j in i]) for i in zip(data.height.values,data.density.values,data.overhang.values)]
#vals = [tuple([j for j in i]) for i in zip(data['height'],data['density'],data['k'],data['overhang'])]
# ...

############################
# GetClimateData # 
inputdir = "/Users/rmichie/Desktop/2003_PRE_CCC/heatsource_9/inputfiles/"
climatefiles = ['Climate_01.csv','Climate_02.csv','Climate_03.csv','Climate_04.csv']
climatedata = pd.DataFrame()

for file in climatefiles:
    newfile = csv_reader(inputdir,file)
    newfile = pd.read_csv(join(inputdir,file.strip()),quotechar='"',quoting=0, index_col='DateTime')
    climatedata = pd.concat([climatedata,newfile],axis=1)    
    
    #newfile = pd.read_csv(join(IniParams["inputdir"],file.strip()),quotechar='"',quoting=0, index_col='DateTime')
    #climatedata = pd.concat([climatedata,newfile],axis=1)
    
data = [tuple(zip(line[0:None:4],line[1:None:4],line[2:None:4],line[3:None:4])) for line in climatedata.values]

x = x
