import csv
import pandas as pd
from bs4 import BeautifulSoup
import re
#import pprint
import urllib
from urllib.request import urlopen

#display list
def display_list(arr):
    return "; ".join(arr)

#get UniProt Protein IDs
input_file = "ProteinIDs_test.csv"
df = pd.read_csv(input_file)

#Get molecular functions
with open('Extracted Data/Protein_data.csv', 'w') as csvfile:
    datawriter = csv.writer(csvfile, delimiter=',')
    header =    [  "Protein ID",
                    "Molecular Function-GO Annot","Molecular Function-Keyword",
                    "Biological processes-GO Annot","Biological processes-Keywords",
                    "Cellular Component-Go Annot", "Cellular Component-Keywords",
                    "Disease-OMIM ID", "Disease-Keywords",
                    "Technical Terms-Keywords"
                ]
    datawriter.writerow(header)

    for i in df['IDs']:

        print(i)
        #arrays for data collection
        molecular_function_go = []
        molecular_function_kw = []
        biological_process_go = []
        biological_process_kw = []
        cellular_component_go = []
        cellular_component_kw = []

        #specify the url
        url = "http://www.uniprot.org/uniprot/" + str(i)

        #Query the website
        page = urlopen(url)

        #Parse the html, store it in Beautiful Soup format
        bsf = BeautifulSoup(page, "lxml")

        #Molecular Function GO Annotation
        ext_data = bsf.find('ul', class_='noNumbering molecular_function')
        for data in ext_data.findAll('li'):
            cells = data.find(lambda tag: tag.name == "a" and (tag.has_attr("onclick")))
            if(not(cells is None)):
                cells = cells.find(text=True).strip()
                molecular_function_go.append(cells)

        #Biological Processes GO Annotation
        ext_data = bsf.find('ul', class_='noNumbering biological_process')
        for data in ext_data.findAll('li'):
            cells = data.find(lambda tag: tag.name == "a" and (tag.has_attr("onclick")))
            if(not(cells is None)):
                cells = cells.find(text=True).strip()
                biological_process_go.append(cells)

        #Cellular Component GO Annotation
        ext_data = bsf.find('ul', class_='noNumbering subcellLocations')
        for data in ext_data.findAll('li'):
            data1 = data.findAll('ul')
            for data2 in data1:
                data3 = data2.findAll('li')
                for data4 in data3:
                    #print(data4.find(text=True))
                    cells = data4.findAll(lambda tag: tag.name == "a" and (tag.has_attr("href")))
                    for cell in cells:
                        c = str(cell)
                        '''
                        if(c.find("locations") != -1):
                            print(cell)
                            #print(cell.find(text=True))
                        '''

        #Keywords - Molecular Function and Biological processes
        ext_data = bsf.find('table', class_='databaseTable')
        ext_data = ext_data.findAll('tr')
        for row in ext_data:
            data = row.findAll('td')
            head = data[0].find(text=True)
            val = data[1].findAll(text=True)
            val = list(filter(lambda x : x != ', ', val))
            if(head=="Molecular function"):
                molecular_function_kw = val
            if(head=="Biological process"):
                biological_process_kw = val

        #write data to CSV file
        datawriter.writerow([i, display_list(molecular_function_go), display_list(molecular_function_kw), display_list(biological_process_go), display_list(biological_process_kw), display_list(cellular_component_go)])
