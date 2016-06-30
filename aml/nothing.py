import sqlite3 as sql

worksheet = '160628_merged'

samples = ['D15-18331', 'D15-21584', 'D15-22373', 'D15-20343', 'D15-25430', 'D03-21521', 'D15-08791', 'D15-08798',
           'D15-04183', 'D15-00899', 'D15-50424', 'D15-45066', 'D14-45300', 'D15-35262', 'D14-33938', 'D13-42537',
           'D14-16565', 'D15-31492', 'D15-41762', 'D14-30832', 'D14-27112', 'D15-02217', 'D15-26810'
           ]

for sample in samples:
    con = sql.connect('/home/cuser/PycharmProjects/django_apps/mypipeline/db.sqlite3')
    curs = con.cursor()
    curs.execute("INSERT INTO Samples (sample, run) VALUES (?,?)", (sample, worksheet))
    con.commit()
