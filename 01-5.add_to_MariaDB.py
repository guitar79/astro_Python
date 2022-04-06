import cv2
import os
import numpy as np
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
db_host = '10.114.0.120'

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

base_drs = ['../CCD_obs_raw/']

fullnames = []
for dirName in base_drs :
    #dirName = "../Aerosol/MODIS Aqua C6.1 - Aerosol 5-Min L2 Swath 3km/2002/185/"
    try :
        fullnames.extend(Python_utilities.getFullnameListOfallFiles("{}".format(dirName)))
    except Exception as err :
        #Python_utilities.write_log(err_log_file, err)
        print(err)
        continue
fullnames = sorted(fullnames)
#########################################

for fullname in fullnames :
    #fullname = fullnames[10]
    if fullname[-4:].lower() == ".fit":
        print("Starting: {}".format(fullname))
        fullname_el = fullname.split("/")
        filename_el = fullname_el[-1].split("/")

        hdul = fits.open(fullname)

        if 'OBJCTRA' in hdul[0].header and \
            'OBJCTDEC' in hdul[0].header and \
            'OBJCTALT' in hdul[0].header and \
            'OBJCTAZ' in hdul[0].header and \
            'OBJCTHA' in hdul[0].header and \
            'EXPTIME' in hdul[0].header :
            try:
                qry = """SELECT `id` FROM `{}`.`{}` WHERE `fullname`= '{}';""".format(db_name, tb_name, fullname)
                check_row = cur.execute(qry)
                print("check_row: {}".format(check_row))

                if check_row == 0:
                    qry = """INSERT INTO `{0}`.`{1}`                                         
                                        (`fullname`, `OBJCTRA`, `OBJCTDEC`, `OBJCTALT`, `OBJCTAZ`, `OBJCTHA`, `EXPTIME`) 
                                        VALUES ('{2}', '{3}', '{4}', '{5}', '{6}', '{7}', '{8}');""".format(db_name, tb_name,
                                                fullname, hdul[0].header['OBJCTRA'], hdul[0].header['OBJCTDEC'],
                                                          hdul[0].header['OBJCTALT'], hdul[0].header['OBJCTAZ'],
                                                          hdul[0].header['OBJCTHA'], hdul[0].header['EXPTIME'])
                    print("qry: {}".format(qry))

                else:
                    qry = """UPDATE `{0}`.`{1}` 
                            SET `OBJCTRA`= '{3}' ,
                            `OBJCTDEC`= '{4}' ,
                            `OBJCTALT`= '{5}' ,
                            `OBJCTAZ`= '{6}' ,
                            `OBJCTHA`= '{7}' ,
                            `EXPTIME`= '{8}'   
                            WHERE `{1}`.`id` = {2};""".format(db_name, tb_name,
                                                fullname, hdul[0].header['OBJCTRA'], hdul[0].header['OBJCTDEC'],
                                                          hdul[0].header['OBJCTALT'], hdul[0].header['OBJCTAZ'],
                                                          hdul[0].header['OBJCTHA'], hdul[0].header['EXPTIME'])
                    print("qry: {}".format(qry))
                cur.execute(qry)
                conn.commit()
            except Exception as err :
                #Python_utilities.write_log(err_log_file, err)
                print("err with {}".format(fullname))
                continue