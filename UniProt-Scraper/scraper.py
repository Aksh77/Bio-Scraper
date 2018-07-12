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
input_file = "UniProt-Scraper/ProteinIDs.csv"
df = pd.read_csv(input_file)
IDs = df['IDs']

#Get molecular functions
with open('UniProt-Scraper/Extracted Data/Protein_data.csv', 'w') as csvfile1, open('UniProt-Scraper/Extracted Data/Comp_bias.csv', 'w') as csvfile2, open('UniProt-Scraper/Extracted Data/Secondary_structures.csv', 'w') as csvfile3, open('UniProt-Scraper/Extracted Data/3D_structures.csv', 'w') as csvfile4:

    #Write header row for Protein_data file
    datawriter1 = csv.writer(csvfile1, delimiter=',')
    header1 =    [  "Protein ID", "Gene ID", "CCDS ID", "BioGrid ID",
                    "Molecular Function-GO Annot","Molecular Function-Keyword",
                    "Biological processes-GO Annot","Biological processes-Keywords",
                    "Cellular Component-Go Annot", "Cellular Component-Keywords",
                    "Disease-OMIM ID", "Disease-Keywords",
                    "Technical Terms-Keywords", "Polymorphism" ]
    datawriter1.writerow(header1)

    #Write header row for Comp_bias file
    datawriter2 = csv.writer(csvfile2, delimiter=',')
    header2 =    [  "Protein ID",
                    "Poly Repeats-Position", "Poly Repeats-Amino Acids", "Poly Repeats-Length",
                    "Rich Domain-Position", "Rich Domain-Amino Acid" ]
    datawriter2.writerow(header2)

    #Write header row for Secondary structure
    datawriter3 = csv.writer(csvfile3, delimiter=',')
    header3 = ["Protein ID", "Turns", "Beta-Strands", "Helixes"]
    datawriter3.writerow(header3)

    #Write header row for 3D structure
    datawriter4 = csv.writer(csvfile4, delimiter=',')
    header4 = ["Protein ID",
                "X-Ray PDB ID", "X-Ray Resolution", "X-Ray Position",
                "NMR PDB ID", "NMR Resolution", "NMR Position"]
    datawriter4.writerow(header4)

#***********************extract data for each protein ID********************************************
    countdown = len(IDs) #Because I am impatient!
    for i in IDs:

        countdown -= 1
        print(countdown, i)
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
        if(not (ext_data is None)):
            data = ext_data.findAll(text=True)
            if "GeneID" in data:
                i = data.index("GeneID")
                gid = data[i+2]

        #get CCDS ID
        ccds = ""
        ext_data = bsf.find('table', class_='databaseTable SEQUENCE')
        if(not (ext_data is None)):
            data = ext_data.findAll(text=True)
            if "CCDS" in data:
                i = data.index("CCDS")
                ccds = data[i+2]

        #get BioGrid ID
        biogrid = ""
        ext_data = bsf.find('table', class_='databaseTable INTERACTION')
        if(not (ext_data is None)):
            data = ext_data.findAll(text=True)
            if "BioGrid" in data:
                i = data.index("BioGrid")
                biogrid = data[i+2]

        #Molecular Function GO Annotation
        molecular_function_go = []
        ext_data = bsf.find('ul', class_='noNumbering molecular_function')
        if(not (ext_data is None)):
            for data in ext_data.findAll('li'):
                cells = data.find(lambda tag: tag.name == "a" and (tag.has_attr("onclick")))
                if(not(cells is None)):
                    cells = cells.find(text=True).strip()
                    molecular_function_go.append(cells)

        #Biological Processes GO Annotation
        biological_process_go = []
        ext_data = bsf.find('ul', class_='noNumbering biological_process')
        if(not (ext_data is None)):
            for data in ext_data.findAll('li'):
                cells = data.find(lambda tag: tag.name == "a" and (tag.has_attr("onclick")))
                if(not(cells is None)):
                    cells = cells.find(text=True).strip()
                    biological_process_go.append(cells)

        #Cellular Component GO Annotation
        cellular_component_go = []
        ext_data = bsf.find('div', id='table-go_annotation')
        if(not (ext_data is None)):
            lt = ext_data.find('ul', class_='noNumbering subcellLocations')
            for li in lt.findAll('li'):
                cell = li.find('h6')
                if(not(cell is None)):
                    cell_loc = cell.find(text=True)
                    cellular_component_go.append(cell_loc)

        #cellular Component Keywords
        cellular_component_kw = []
        ext_data = bsf.find('div', class_='section ', id="subcellular_location")
        if(not (ext_data is None)):
            header = ext_data.find('h4')
            if(not (header is None)):
                head_data = header.findAll(text=True)
                if 'Keywords - Cellular component' in head_data:
                    data = header.next_sibling
                    vals = data.findAll('a')
                    val = [v.find(text=True) for v in vals]
                    v = list(filter(lambda x : x != ', ', vals))
                    val = []
                    for kw in v:
                        kws = str(kw)
                        if kws[:18]=='<a href="/keywords':
                            val.append(kw.find(text=True))
                    cellular_component_kw = val

        #Keywords - Molecular Function and Biological processes
        molecular_function_kw = []
        biological_process_kw = []
        ext_data = bsf.find('table', class_='databaseTable')
        if(not (ext_data is None)):
            ext_data = ext_data.findAll('tr')
            for row in ext_data:
                data = row.findAll('td')
                head = data[0].find(text=True)
                vals = data[1].findAll('a')
                val = [v.find(text=True) for v in vals]
                v = list(filter(lambda x : x != ', ', vals))
                val = []
                for kw in v:
                    kws = str(kw)
                    if kws[:18]=='<a href="/keywords':
                        val.append(kw.find(text=True))
                if(head=="Molecular function"):
                    molecular_function_kw = val
                if(head=="Biological process"):
                    biological_process_kw = val

        #Disease OMIM ID
        disease_omim_id = []
        ids = []
        ext_data = bsf.findAll('div', class_='diseaseAnnotation')
        if(not (ext_data is None)):
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
        if(not (ext_data is None)):
            heads = ext_data.findAll('h4')
            for head in heads:
                data = head.findAll(text=True)
                if 'Keywords - Disease' in data:
                    j = data.index('Keywords - Disease')
                    val = data[j].parent.parent
                    cells = val.next_sibling
                    vals = cells.findAll('a')
                    val = [v.find(text=True) for v in vals]
                    v = list(filter(lambda x : x != ', ', vals))
                    val = []
                    for kw in v:
                        kws = str(kw)
                        if kws[:18]=='<a href="/keywords':
                            val.append(kw.find(text=True))
                    disease_kw = val
                    break

        #Technical Terms - Keywords
        tech_term_kw = []
        ext_data = bsf.find('div', class_='section', id='miscellaneous')
        if(not (ext_data is None)):
            heads = ext_data.findAll('h4')
            for head in heads:
                data = head.findAll(text=True)
                if 'Keywords - Technical term' in data:
                    j = data.index('Keywords - Technical term')
                    val = data[j].parent.parent
                    cells = val.next_sibling
                    vals = cells.findAll('a')
                    val = [v.find(text=True) for v in vals]
                    v = list(filter(lambda x : x != ', ', vals))
                    val = []
                    for kw in v:
                        kws = str(kw)
                        if kws[:18]=='<a href="/keywords':
                            val.append(kw.find(text=True))
                    tech_term_kw = val
                    break

        #Polymorphism
        polymorphism = ""
        ext_data = bsf.find('div', class_='section', id='sequences')
        if(not (ext_data is None)):
            heads = ext_data.findAll('h4')
            for head in heads:
                data = head.findAll(text=True)
                if 'Polymorphism' in data:
                    j = data.index('Polymorphism')
                    val = data[j].parent.parent
                    polymorphism = val.next_sibling.find(text=True)
                    break

        #get compositional bias data
        poly_pos = []
        poly_desc = []
        poly_len = []
        rich_pos = []
        rich_desc = []
        ext_data = bsf.find('table', class_='featureTable', id="Compositional_bias_section")
        if(not(ext_data is None)):
            for row in ext_data.findAll('tr'):
                num = row.findAll('td', class_='numeric')
                if(len(num)>0):
                    pos = num[0].find(text=True)
                    pos = '-'.join(pos.split('\xa0–\xa0'))
                    length = num[1].find(text=True)
                    desc = row.findAll('td', class_='featdescription')
                    comp_desc = desc[0].find(text=True)
                    if(comp_desc[:4])=='Poly':
                        poly_pos.append(pos)
                        poly_desc.append(comp_desc[5:])
                        poly_len.append(length)
                    else:
                        rich_pos.append(pos)
                        rich_desc.append(comp_desc[:-5])
        polypos = '\n'.join(poly_pos)
        polydesc = '\n'.join(poly_desc)
        polylen = '\n'.join(poly_len)
        richpos = '\n'.join(rich_pos)
        richdesc = '\n'.join(rich_desc)

        #get Secondary structures data
        helix = []
        beta_strand = []
        turn = []

        ext_data = bsf.find('table', class_='featureTable', id="secstructure_section")
        if(not(ext_data is None)):
            for row in ext_data.findAll('tr'):
                cols = row.findAll('td')
                if len(cols)>0:
                    desc = cols[0].findAll(text=True)
                    desc = desc[1]
                    pos = cols[1].find(text=True)
                    pos = re.findall(r'\d+', pos)
                    pos = '-'.join(pos)
                    if desc == 'Helix':
                        helix.append(pos)
                    elif desc == 'Beta strand':
                        beta_strand.append(pos)
                    elif desc == 'Turn':
                        turn.append(pos)
        helixes = '\n'.join(helix)
        beta_strands = '\n'.join(beta_strand)
        turns = '\n'.join(turn)

        #get 3D structures data
        x_id = []
        x_res = []
        x_pos = []
        n_id = []
        n_res = []
        n_pos = []
        ext_data = bsf.find('table', class_='databaseTable STRUCTURE')
        if(not (ext_data is None)):
            for row in ext_data.findAll('tr'):
                data = row.findAll(text=True)
                if 'X-ray' in data:
                    x_id.append(data[0])
                    x_res.append(data[2])
                    x_pos.append(data[4])
                elif 'NMR' in data:
                    n_id.append(data[0])
                    n_res.append(data[2])
                    n_pos.append(data[4])
        xid = '\n'.join(x_id)
        xres = '\n'.join(x_res)
        xpos = '\n'.join(x_pos)
        nid = '\n'.join(n_id)
        nres = '\n'.join(n_res)
        npos = '\n'.join(n_pos)

#***************************Write extracted data to CSV files**********************************
        #write data to Protein_data CSV file
        datawriter1.writerow([pid, gid, ccds, biogrid,
                            display_list(molecular_function_go), display_list(molecular_function_kw),
                            display_list(biological_process_go), display_list(biological_process_kw),
                            display_list(cellular_component_go), display_list(cellular_component_kw),
                            display_list(disease_omim_id), display_list(disease_kw),
                            display_list(tech_term_kw), polymorphism])

        #write data to Comp_bias CSV file
        datawriter2.writerow([pid, polypos, polydesc, polylen, richpos, richdesc])

        #write data to Secondary_structures CSV file
        datawriter3.writerow([pid, turns, beta_strands, helixes])

        #write data to 3D structures CSV file
        datawriter4.writerow([pid, xid, xres, xpos, nid, nres, npos])
