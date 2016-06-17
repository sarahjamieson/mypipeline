import sqlite3 as sql

con = sql.connect('/home/cuser/PycharmProjects/django_apps/mysitedev/mypipeline.db.sqlite3')
curs = con.cursor()

curs.execute("CREATE TABLE IF NOT EXISTS Results(result_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, run TEXT, "
             "sample TEXT, caller TEXT, chrom TEXT, pos INTEGER, ref TEXT, alt TEXT, chr2 TEXT, end INTEGER, "
             "region TEXT, sv_type TEXT, size INTEGER, gt TEXT, depth INTEGER, alleles TEXT, ab REAL, gene TEXT, "
             "func TEXT, exonic_func TEXT)")
curs.execute("CREATE TABLE IF NOT EXISTS Runs(run_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, run TEXT)")
curs.execute("CREATE TABLE IF NOT EXISTS Samples(sample_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, sample TEXT, "
             "run TEXT")


# along with adding to excel spreadsheet
run = ''
if curs.execute("SELECT * FROM Runs WHERE run LIKE '%s'" % run) is None:
    curs.execute("INSERT INTO Runs (run) VALUES (?)", (run,))
else:
    pass
sample = ''
if curs.execute("SELECT * FROM Samples WHERE sample LIKE '%s'" % sample) is None:
    curs.execute("INSERT INTO Samples (sample) VALUES (?)", (sample,))
else:
    pass

curs.execute("INSERT INTO Results (run, sample, caller, chrom, pos, ref, alt, chr2, end, region, sv_type, size, gt, "
             "depth, alleles, ab, gene, func, exonic_func) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
             (run, sample, caller, chrom, pos, ref, alt, chr2, end, region, sv_type, size, gt, depth, alleles, ab,
              gene, func, exonic_func))
con.commit()
