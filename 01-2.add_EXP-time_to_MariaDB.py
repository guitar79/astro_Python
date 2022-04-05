import cv2
import os
import pandas as pd 
import numpy as np
import exifread
import pymysql
from urllib import parse
import Python_utilities
import FB107_utilities

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

tb_name = 'SAVE_file_info'

conn = pymysql.connect(host=db_host, port=3306,
                       user=db_user, password=db_pass, db=db_name,
                       charset='utf8mb4', cursorclass=pymysql.cursors.DictCursor)

cur = conn.cursor()

qry = "SHOW TABLES FROM `{}` LIKE '{}'".format(db_name, tb_name)

table_exist = cur.execute(qry)

if table_exist == 0:        
    qry = """CREATE TABLE IF NOT EXISTS `{}`.`{}` (
        `id` INT NOT NULL AUTO_INCREMENT ,
        `fullname` VARCHAR(16384) default NULL ,
        `brightness-mean` VARCHAR(16) default NULL ,
        `brightness-std` VARCHAR(16) default NULL ,
        `Line` VARCHAR(5) default NULL ,
    	`EXIF` TEXT default NULL ,
    	`EXP-time` VARCHAR(10) NULL DEFAULT NULL AFTER `EXIF`;
        `Update-DT` TIMESTAMP on update CURRENT_TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP , 
        PRIMARY KEY (`id`)) ENGINE = InnoDB;""".format(db_name, tb_name)
    cur.execute(qry)
    conn.commit()
else :
    qry = "SELECT COLUMN_NAME FROM INFORMATION_SCHEMA.COLUMNS WHERE TABLE_SCHEMA = '{}' AND TABLE_NAME = '{}';".format(db_name, tb_name)
    cur.execute(qry)
    colunm_dic = cur.fetchall()
    print(colunm_dic)

#########################################

qry = "SELECT `fullname` FROM `{}`.`{}` WHERE `EXP-time` IS NULL;".format(db_name, tb_name)
cur.execute(qry)
fullnames = cur.fetchall()

for fullname in fullnames :
    print(fullname['fullname'])
    with open(fullname['fullname'], 'rb') as f:
        tags = exifread.process_file(f)
        print(tags["EXIF ExposureTime"])
        qry = """UPDATE `{0}`.`{1}` 
                SET `EXP-time` = '{2}'   
                WHERE `{1}`.`fullname` = '{3}';""".format(db_name, tb_name,
                                                  tags["EXIF ExposureTime"], fullname['fullname'])
        print("q3_update: {}".format(qry))
        cur.execute(qry)
    conn.commit()