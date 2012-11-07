#!/usr/bin/env python
###############################################################################
#
#    CSVtoSqlite.py
#    
#    convert our Casy CSVs to sqlite for further parsing
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
from pprint import pprint
import sqlite3 as sqlite
import os
import csv
import string
from Bio import SeqIO
import glob

###############################################################################
# CODE HERE
###############################################################################

def doWork( options ):
    # connect to the DB
    dbfname = options.sqlite_file
    if os.path.isfile(dbfname):
        os.remove(dbfname)
    db = sqlite.connect(dbfname)
    
    # create the table and columns we'll be needing
    cur = db.cursor()
    cur.execute('CREATE TABLE '+options.table_name+' (gbid TEXT PRIMARY KEY, cazy_family TEXT, domain TEXT, prot_name TEXT, prot_desc TEXT, organism TEXT, tax TEXT, ec TEXT, seq TEXT, comment TEXT)')
    db.commit()
    
    # read in the CSV files one by one and insert their elements into the table
    dirList=glob.glob(options.procpath+'/*.csv')
    for csv_fn in dirList:
        csv_name_fields = os.path.basename(csv_fn).split('.')
        cazy_family = csv_name_fields[0] 
        print "Loading data from", csv_fn
        
        # open the input and output files
        csv_fh = csv.reader(open(csv_fn, 'rb'), delimiter='\t', quotechar="'")
        # skip header
        csv_fh.next()
        for row in csv_fh:
            # check to see if this mofo is unique
            cur.execute('SELECT COUNT(*) FROM '+options.table_name+' WHERE gbid="'+row[4]+'"')
            (number_of_rows,)=cur.fetchone()
            if(number_of_rows == 0):
                # unique!
                # check to see if the GB file is in the path
                gbfile = os.path.join(options.GB_path, row[4]+".gb")
                if os.path.isfile(gbfile):
                    prot_name = row[1]
                    prot_desc = ""
                    organism = row[2]
                    tax = ""
                    ec_num = row[3]
                    seq = ""
                    comment = ""
                    # only one record per file!
                    for seq_record in SeqIO.parse(gbfile, "genbank"):
                        seq = str(seq_record.seq)
                        annotations = seq_record.annotations
                        if("" == organism):
                            if('organism' in annotations):
                                organism = annotations['organism']
                        if('comment' in annotations):
                            comment = annotations['comment']
                        if('taxonomy' in annotations):    
                            tax = string.join(annotations['taxonomy'], ",")
                    
                    # extract the description from the features
                    for feat in seq_record.features:
                        if(feat.type == 'Protein'):
                            if('product' in feat.qualifiers):
                                prot_desc = string.join(feat.qualifiers['product'], ",")

                    # strip quotes!
                    prot_desc = stripQuotes(prot_desc)
                    prot_name = stripQuotes(prot_name)
                    tax = stripQuotes(tax)
                    comment = stripQuotes(comment)
                    organism = stripQuotes(organism)
                    
                    try:
                        cur.execute('INSERT INTO '+options.table_name+' (gbid, cazy_family, domain, prot_name, prot_desc, organism, tax, ec, seq, comment) VALUES("'+row[4]+'", "'+cazy_family+'", "'+row[0]+'", "'+prot_name+'", "'+prot_desc+'", "'+organism+'", "'+tax+'", "'+ec_num+'", "'+seq+'", "'+comment+'")')
                    except Exception:
                        print [row[4],cazy_family,row[0],prot_name,prot_desc,organism,tax,ec_num,seq,comment]
                        raise
                    db.commit()
        
    return 0

def stripQuotes(inStr):
    return string.replace(string.replace(inStr, "'", ""), "\"", "")

###############################################################################
# TEMPLATE SUBS
###############################################################################
#
# Entry point, parse command line args and call out to doWork
#
if __name__ == '__main__':
    # intialise the options parser
    parser = argparse.ArgumentParser(description='Convert notrmalised Casy CSVs to an sqllite DB')
    
    # add positional options here:
    parser.add_argument('procpath', help="Directory containing processed csv files to parse")
    parser.add_argument('GB_path', help="SQLDB name")
    parser.add_argument('table_name', help="Name of this table in the SQL lite DB")
    parser.add_argument('sqlite_file', help="SQLDB name")
    
    # add optional options here
    #parser.add_argument('-m', '--max_rows', default=0, help="Load only this many rows of data [default: load all rows]")
    #parser.add_argument('-p', '--plot', action="store_true", default=False, help="Plot the PCA")
    
    # get and check options
    args = parser.parse_args()

    # 
    # do what we came here to do
    #
    doWork(args)
