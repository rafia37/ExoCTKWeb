import os
import sqlite3
import datetime

def create_db(dbpath, schema, overwrite=True):
    """
    Create a new database at the given dbpath

    Parameters
    ----------
    dbpath: str
        The full path for the new database, including the filename and .db file extension.
    schema: str
        The path to the .sql schema for the database
    overwrite: bool
        Overwrite dbpath if it already exists
    """
    if dbpath.endswith('.db'):
        
        if os.path.isfile(dbpath) and overwrite:
            os.system('rm {}'.format(dbpath))
            
        os.system("cat {} | sqlite3 {}".format(schema, dbpath))
            
        if os.path.isfile(dbpath):
            print("ExoCTK database created at {}".format(dbpath))
        
    else:
        print("Please provide a path and file name with a .db file extension, e.g. /Users/<username>/Desktop/test.db")

def load_db(dbpath):
    if os.path.isfile(dbpath):
    
        con = sqlite3.connect(dbpath, isolation_level=None, detect_types=sqlite3.PARSE_DECLTYPES, check_same_thread=False)
        cur = con.cursor()
        
        return cur
        
    else:
        print("Sorry, could not find the file '{}'".format(dbpath))

def log_form_input(form_dict, table, database):
    """
    A function to store the form inputs of any page GET requests in a database
    """
    # Convert hyphens to underscores and leading numerics to letters for db column names
    inpt = {k.replace('-','_').replace('3','three').replace('4','four'):v[0] for k,v in dict(form_dict).items()}
    
    # Add a timestamp
    inpt['date'] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    
    # Insert the form valsues
    qmarks = ', '.join('?' * len(inpt))
    qry = "Insert Into {} ({}) Values ({})".format(table, ', '.join(inpt.keys()), qmarks)
    database.execute(qry, list(inpt.values()))