# -- Import libraries -- # 

import urllib.parse
import urllib.request   
import gzip                  

# -- Functions -- #

def map_db_id_uniprot(query, from_db = "ID", to_db = "GENENAME"):
    """
    Uses the UniProt API to map IDs from different databases. For more 
    information: https://www.uniprot.org/help/api_idmapping

    Parameters
    ----------
    from_db: str
           query's original database id 
    to_db: str
           query's desired database id 
    query: str 
           list of gene/protein names as a string (comma-separated) 
    
    
    Returns
    -------
    response.split("\t")[-1].strip(): str
            Gene/protein name in the desired database id
    """

    url = 'https://www.uniprot.org/uploadlists/'

    params = {
    'from': from_db,
    'to': to_db,
    'format': 'tab',
    'query': query
    }

    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)

    with urllib.request.urlopen(req) as f:
        response = f.read().decode('utf-8')
    
    return response.split("\t")[-1].strip()
       

def map_ac_name_proteome(proteome_dir):
       """
       Generates a dictionary with the proteomes correspondance
       between AC and gene name

       Parameters
       ----------
       proteome_dir: str
           Path to the folder where the compressed proteome FASTA file is located
    
    
       Returns
       -------
       ac_gene_name: dict
            Dictionary containing each protein's accession number as key and the gene's name as value
       """

       with gzip.open(proteome_dir, 'rt') as f:         # important to set read mode as 'rt' because we are reading a compressed text file
        
              data = f.readlines()
              ac_gene_name_l = []  # dicts cannot be created as dict[key] = value
                                   # if the key is a int, so we are creating it from
                                   # a list

              for line in data:

                     if line[0] == '>': # Description lines start with “>” character
                            ac = line.split()[0][1:].split("|")[1] # The first column of the description line is the sequence id
                                                                      # the second element of such description is the AC 
                            gene_name = line.split()[0][1:].split("|")[2].split("_")[0]    # the third is the gene name
                            ac_gene_name_l.append((ac, gene_name))

       
       return dict(ac_gene_name_l)
       


                     