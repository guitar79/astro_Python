import os
from astropy.io import fits
import Python_utilities

#########################################
log_dir = "logs/"
log_file = "{}{}.log".format(log_dir, os.path.basename(__file__)[:-3])
err_log_file = "{}{}_err.log".format(log_dir, os.path.basename(__file__)[:-3])
print ("log_file: {}".format(log_file))
print ("err_log_file: {}".format(err_log_file))
#########################################

#########################################
# mariaDB info
#########################################
import pymysql
db_host = '192.168.0.20'
#db_host = '10.114.0.120'

db_user = 'modis'
db_pass = 'Modis12345!'
db_name = 'MODIS_Aerosol'

db_user = 'FB107'
db_pass = 'Gses12345!'
db_name = 'FB107'

tb_name = 'FITS_file_info'

conn = pymysql.connect(host=db_host, port=3306,
                       user=db_user, password=db_pass, db=db_name,
                       charset='utf8mb4', cursorclass=pymysql.cursors.DictCursor)

cur = conn.cursor()

qry = "SHOW TABLES FROM `{}` LIKE '{}'".format(db_name, tb_name)

table_exist = cur.execute(qry)

if table_exist == 0:
    qry = """CREATE TABLE IF NOT EXISTS `{}`.`{}` (
        `id` INT NOT NULL AUTO_INCREMENT ,
        `fullname` VARCHAR(16384) NULL default NULL ,
        `date-time` VARCHAR(14) NULL default NULL ,
        `Frame_TYP` VARCHAR(12) NULL DEFAULT NULL ,
        `OBJCTRA` VARCHAR(10) NULL DEFAULT NULL ,      
        `OBJCTDEC` VARCHAR(10) NULL DEFAULT NULL ,          
        `OBJCTALT` VARCHAR(10) NULL DEFAULT NULL ,             
        `OBJCTAZ` VARCHAR(10) NULL DEFAULT NULL ,              
        `OBJCTHA` VARCHAR(10) NULL DEFAULT NULL ,   
    	`EXPTIME` VARCHAR(20) NULL DEFAULT NULL ,
        `Update-DT` TIMESTAMP on update CURRENT_TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP , 
        PRIMARY KEY (`id`)) ENGINE = InnoDB;""".format(db_name, tb_name)
    cur.execute(qry)
    conn.commit()
else :
    qry = "SELECT COLUMN_NAME FROM INFORMATION_SCHEMA.COLUMNS WHERE TABLE_SCHEMA = '{}' AND TABLE_NAME = '{}';".format(db_name, tb_name)
    cur.execute(qry)
    column_dic = cur.fetchall()
    print(column_dic)

#########################################

qry = "SELECT `fullname` FROM `{}`.`{}` WHERE `date-time` IS NULL;".format(db_name, tb_name)
cur.execute(qry)
fullnames = cur.fetchall()

for fullname in fullnames :
    #fullname = fullnames[0]
    print(fullname['fullname'])
    fullname_el = fullname['fullname'].split("/")
    filename_el = fullname_el[-1].split("_")
    d_time = filename_el[3]
    d_time = d_time.replace("-", "")

    qry = """UPDATE `{0}`.`{1}` 
            SET `date-time`= '{3}'    
            WHERE `{1}`.`fullname` = '{2}';""".format(db_name, tb_name,
                                              fullname['fullname'], d_time)
    print("qry: {}".format(qry))
    cur.execute(qry)
    conn.commit()

