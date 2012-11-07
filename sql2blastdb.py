#!/usr/bin/env python
###############################################################################
#
#    CSVtoSqlite.py
#    
#    reach into the sqlite db and make a blast db
#
#    Copyright (C) 2012 Michael Imelfort
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

import argparse
import sys
import sqlite3 as sqlite
import os
import string
from argparse import RawTextHelpFormatter
import time

###############################################################################
# CODE HERE
###############################################################################

def doWork( options ):
    # connect to the DB
    db = sqlite.connect(options.sqlite_file)

    # get the single table name
    fn_fields = os.path.basename(options.sqlite_file).split('.')
    table_name = fn_fields[0] 
    cur = db.cursor()
    
    # work out the search clause
    where_clause = ""
    if(options.query):
        where_clause = " WHERE "+ options.query
    query = 'SELECT gbid,prot_desc,seq FROM '+table_name+where_clause

    # do it!
    cur.execute(query)
    db.commit()
    results=cur.fetchall()
    
    # make the fasta file
    FILE = open("tmp.fasta","w")
    for result in results:
        FILE.writelines(">" + result[0] + "_" + result[1] + "\n" + result[2] + "\n")
#        FILE.writelines(">" + result[0] + "\n" + result[1] + "\n")
    FILE.close()
    
    blast_cmd = "makeblastdb -dbtype prot -in tmp.fasta -title "+options.blast_name+" -parse_seqids -out "+options.blast_name
    os.system(blast_cmd)
    
    # clean up
    os.system("rm tmp.fasta")
    return 0

###############################################################################
# TEMPLATE SUBS
###############################################################################
#
# Entry point, parse command line args and call out to doWork
#
if __name__ == '__main__':
    # intialise the options parser
    desc = "Reach into a cazy sqlite db and make a blast db"
    
    epi = """
    Scraped DBs are:
    Glycoside Hydrolases (GH.sqlite) : hydrolysis and/or rearrangement of glycosidic bonds (see CAZypedia definition)
    GlycosylTransferases (GT.sqlite) : formation of glycosidic bonds (see definition)
    Polysaccharide Lyases (PL.sqlite) : non-hydrolytic cleavage of glycosidic bonds
    Carbohydrate Esterases (CE.sqlite) : hydrolysis of carbohydrate esters
    Carbohydrate-Binding Modules (CBM.sqlite) : adhesion to carbohydrates

    Table columns are:
    
    gbid                   Genbank id - contains version (EX: AAD25394.1)
    cazy_family            Cazy sub family (EX GH10)
    domain                 Bacteria,etc...
    prot_name              Code name for the protein (EX: pelA)
    prot_desc              Longer description of the protein (EX: pectate lyase)
    organism               Full organism name   (EX: Azospirillum irakense)
    tax                    Full taxonomy NOTE: May be mixed with organism field 
                (EX: Bacteria,Proteobacteria,Alphaproteobacteria,Rhodospirillales,Rhodospirillaceae,Azospirillum)
    ec                     KEGG EC number
    seq                    Protein sequence
    comment                Comment in genbank file
                (EX: PREDICTED REFSEQ: This record has not been reviewed and the
                function is unknown. The reference sequence was derived from
                EEF92833. Method: conceptual translation.)

    Use SQL like queries to select sequences
    
    Example:
    
    ./sql2blastdb.py PL.sqlite PL_marine.db -q "comment LIKE '%marine%'"
    
    Creates a blastdb called PL_marine.db containing all sequences in the PL database with a comment containing the word "marine"
    
    ./sql2blastdb.py PL.sqlite PL_marine.db -q "tax LIKE '%Rhodospirillales%' AND comment LIKE '%marine%'"
    
    Creates a blastdb  PL_marine.db containing all sequences in the PL database with a comment containing the word "marine"
    And a taxonomy string containing "Rhodospirillales"
    
    Enjoi!
    
    """
    
    parser = argparse.ArgumentParser(description=desc, epilog=epi, formatter_class=RawTextHelpFormatter)
    
    # add positional options here:
    parser.add_argument('sqlite_file', help="SQLDB name")
    parser.add_argument('blast_name', help="Name of the blastdb to create")
    
    # add optional options here
    parser.add_argument('-q', '--query', help="SQL query to parse")
    
    # get and check options
    args = parser.parse_args()

    # 
    # do what we came here to do
    #
    doWork(args)
