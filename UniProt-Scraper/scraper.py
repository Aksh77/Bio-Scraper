import csv
import pandas as pd
from bs4 import BeautifulSoup
import re
import pprint
import urllib
from urllib.request import urlopen

#display list
def display_list(arr):
    return ";\n".join(arr)

#get data; handle error
def fetchdata(url):
    try:
        return urlopen(url)
    except:
        fetchdata(url)

#get UniProt Protein IDs
input_file = "ProteinIDs_test.csv"
df = pd.read_csv(input_file)

#Get molecular functions
with open('Extracted Data/Protein_data.csv', 'w') as csvfile:
    datawriter = csv.writer(csvfile, delimiter=',')
    header =    [  "Protein ID", "Gene ID", "CCDS ID", "BioGrid ID",
                    "Molecular Function-GO Annot","Molecular Function-Keyword",
                    "Biological processes-GO Annot","Biological processes-Keywords",
                    "Cellular Component-Go Annot", "Cellular Component-Keywords",
                    "Disease-OMIM ID", "Disease-Keywords",
                    "Technical Terms-Keywords", "Polymorphism"
                ]
    datawriter.writerow(header)

    for i in df['IDs']:
        print(i)
        pid = i

        #specify the url
        url = "http://www.uniprot.org/uniprot/" + str(i)
        #Query the website
        page = fetchdata(url)
        #Parse the html, store it in Beautiful Soup format
        bsf = BeautifulSoup(page, "lxml")

        #get Gene ID
        gid = ""
        ext_data = bsf.find('table', class_='databaseTable GENOME')
        data = ext_data.findAll(text=True)
        if "GeneID" in data:
            i = data.index("GeneID")
            gid = data[i+2]

        #get CCDS ID
        ccds = ""
        ext_data = bsf.find('table', class_='databaseTable SEQUENCE')
        data = ext_data.findAll(text=True)
        if "CCDS" in data:
            i = data.index("CCDS")
            ccds = data[i+2]

        #get BioGrid ID
        biogrid = ""
        ext_data = bsf.find('table', class_='databaseTable INTERACTION')
        data = ext_data.findAll(text=True)
        if "BioGrid" in data:
            i = data.index("BioGrid")
            biogrid = data[i+2]

        #Molecular Function GO Annotation
        molecular_function_go = []
        ext_data = bsf.find('ul', class_='noNumbering molecular_function')
        for data in ext_data.findAll('li'):
            cells = data.find(lambda tag: tag.name == "a" and (tag.has_attr("onclick")))
            if(not(cells is None)):
                cells = cells.find(text=True).strip()
                molecular_function_go.append(cells)

        #Biological Processes GO Annotation
        biological_process_go = []
        ext_data = bsf.find('ul', class_='noNumbering biological_process')
        for data in ext_data.findAll('li'):
            cells = data.find(lambda tag: tag.name == "a" and (tag.has_attr("onclick")))
            if(not(cells is None)):
                cells = cells.find(text=True).strip()
                biological_process_go.append(cells)

        #Cellular Component GO Annotation
        cellular_component_go = []
        '''
        ext_data = bsf.find('ul', class_='noNumbering subcellLocations')
        for data in ext_data.findAll('li'):
            head = data.find('h6')
            if(not(head is None)):
                data1 = head.next_sibling
                val = data1.findAll('li')
                if(not(val is None)):
                    for cell in val:
                        val = cell.find(text=True)
                        if val[:4]!="Ref.":
                            cellular_component_go.append(val)
        '''

        #cellular Component Keywords
        cellular_component_kw = []
        ext_data = bsf.find('div', class_='section ', id="subcellular_location")
        header = ext_data.find('h4')
        head_data = header.findAll(text=True)
        if 'Keywords - Cellular component' in head_data:
            data = header.next_sibling
            val = data.findAll(text=True)
            val = list(filter(lambda x : x != ', ', val))
            cellular_component_kw = val

        #Keywords - Molecular Function and Biological processes
        molecular_function_kw = []
        biological_process_kw = []
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

        #Disease OMIM ID
        disease_omim_id = []
        ids = []
        ext_data = bsf.findAll('div', class_='diseaseAnnotation')
        for data in ext_data:
            val = data.findAll('a')
            for data1 in val:
                data2 = data.findAll(text=True)
                for j in data2:
                    if j[:8]=="See also":
                        ids.append(j[14:])
        ids = set(ids)
        disease_omim_id = list(ids)

        #Disease Keywords
        disease_kw = []
        ext_data = bsf.find('div', class_='section', id='pathology_and_biotech')
        heads = ext_data.findAll('h4')
        for head in heads:
            data = head.findAll(text=True)
            if 'Keywords - Disease' in data:
                j = data.index('Keywords - Disease')
                val = data[j].parent.parent
                cells = val.next_sibling.findAll(text=True)
                val = list(filter(lambda x : x != ', ', cells))
                disease_kw = val
                break

        #write data to CSV file
        datawriter.writerow([pid, gid, ccds, biogrid,
                            display_list(molecular_function_go), display_list(molecular_function_kw),
                            display_list(biological_process_go), display_list(biological_process_kw),
                            display_list(cellular_component_go), display_list(cellular_component_kw),
                            display_list(disease_omim_id), display_list(disease_kw),
                            ])
