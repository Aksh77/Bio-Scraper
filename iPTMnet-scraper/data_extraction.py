import csv
import pandas as pd
from bs4 import BeautifulSoup
import urllib
from urllib.request import urlopen

#bypass the urllib and “SSL: CERTIFICATE_VERIFY_FAILED” Error
import requests
requests.packages.urllib3.disable_warnings()
import ssl
try:
    _create_unverified_https_context = ssl._create_unverified_context
except AttributeError:
    pass
else:
    ssl._create_default_https_context = _create_unverified_https_context

#display sites
def display_sites(arr):
    l = len(arr)
    arr.sort()
    data = ""
    for i in range(l-1):
        data = data + arr[i] + " "
    data = data + arr[l-1]
    return data

#get iPTMnet IDs for querying
input_file = "ProteinIDs_test.csv"
df = pd.read_csv(input_file)
data = {}

with open('Extracted Data/PTMdata.csv', 'w') as csvfile1, open('Extracted Data/Enzymedata.csv', 'w') as csvfile2:
    datawriter1 = csv.writer(csvfile1, delimiter=',')
    datawriter1.writerow(["Protein ID", "PTM Type", "PTM Site"])
    datawriter2 = csv.writer(csvfile2, delimiter=',')
    datawriter2.writerow(["Protein ID", "PTM Type", "PTM Site", "Enzyme"])
    for i in df['IDs']:
        #dictionary to store data for each iPTMnet ID
        data[i] = {}
        flag = 0
        flag_e = 0

        #specify the url
        url = "https://research.bioinformatics.udel.edu/iptmnet/entry/" + str(i)

        #Query the website
        page = urlopen(url)

        #Parse the html, store it in Beautiful Soup format
        bsf = BeautifulSoup(page, "lxml")

        #find the iptm data table
        table = bsf.find('table', class_='iptm-entry-table')
        tdata = table.find('tbody')

        #get PTM sites and PTM types
        for row in tdata.findAll("tr"):
            cells = row.findAll('td')
            PTMsite = cells[1].find(text=True).strip()
            PTMtype = cells[2].find(text=True).strip()
            enzyme = cells[3].findAll(text=True)
            for j in range(len(enzyme)):
                enzyme[j] = enzyme[j].replace(',','')
                enzyme[j] = enzyme[j].strip()
            enzyme = ' '.join(enzyme)
            #store in dictionary
            if PTMtype in data[i].keys():
                data[i][PTMtype].append(PTMsite)
            else:
                data[i][PTMtype] = [PTMsite]

            #write enzyme data to CSV file
            if flag_e == 0:
                datawriter2.writerow([i, PTMsite, PTMtype, enzyme])
                flag_e = 1
            else:
                datawriter2.writerow(["", PTMsite, PTMtype, enzyme])

        #write PTM data to CSV file
        for k, v in data[i].items():
            if flag==0:
                datawriter1.writerow([i, k, display_sites(v)])
                flag = 1
            else:
                datawriter1.writerow(["", k, display_sites(v)])
        print(i)
