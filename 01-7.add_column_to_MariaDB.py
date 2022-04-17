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
        `OBJECT` VARCHAR(20) NULL DEFAULT NULL ,
        `IMAGETYP` VARCHAR(10) NULL DEFAULT NULL ,
        `DATE-OBS` VARCHAR(28) NULL default NULL ,
        `FOCALLEN` VARCHAR(10) NULL DEFAULT NULL ,
        `XPIXSZ` VARCHAR(10) NULL DEFAULT NULL ,
        `PIXSCALE` VARCHAR(10) NULL DEFAULT NULL, 
        `OBJCTRA` VARCHAR(10) NULL DEFAULT NULL ,      
        `OBJCTDEC` VARCHAR(10) NULL DEFAULT NULL ,          
        `OBJCTALT` VARCHAR(10) NULL DEFAULT NULL ,             
        `OBJCTAZ` VARCHAR(10) NULL DEFAULT NULL ,              
        `OBJCTHA` VARCHAR(10) NULL DEFAULT NULL ,  
        `CTYPE1` VARCHAR(16) NULL DEFAULT NULL ,  
        `CTYPE2` VARCHAR(16) NULL DEFAULT NULL ,  
        `EQUINOX` VARCHAR(16) NULL DEFAULT NULL ,  
        `LONPOLE` VARCHAR(10) NULL DEFAULT NULL ,  
        `LATPOLE` VARCHAR(10) NULL DEFAULT NULL ,  
        `CRVAL1` VARCHAR(16) NULL DEFAULT NULL ,  
        `CRVAL2` VARCHAR(16) NULL DEFAULT NULL ,  
        `CRPIX1` VARCHAR(16) NULL DEFAULT NULL ,  
        `CRPIX2` VARCHAR(16) NULL DEFAULT NULL ,  
        `CUNIT1` VARCHAR(10) NULL DEFAULT NULL ,                                        
        `CUNIT2` VARCHAR(10) NULL DEFAULT NULL ,  
        `CD1_1` VARCHAR(16) NULL DEFAULT NULL ,  
        `CD1_2` VARCHAR(16) NULL DEFAULT NULL ,  
        `CD2_1` VARCHAR(16) NULL DEFAULT NULL ,  
        `CD2_2` VARCHAR(16) NULL DEFAULT NULL ,  
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

ins_column = "OBJECT"
ins_column = "FOCALLEN"
ins_column = "XPIXSZ"
ins_column = "PIXSCALE"

ins_columns = ['IMAGETYP', 'OBJECT', 'FOCALLEN', 'XPIXSZ', 'PIXSCALE', 
                'OBJCTRA', 'OBJCTDEC', 'OBJCTALT', 'OBJCTAZ', 'OBJCTHA', 
                'CTYPE1', 'CTYPE2', 'EQUINOX', 'LONPOLE', 'LATPOLE', 
                'CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2', 'CUNIT1', 
                'CUNIT2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']
#########################################

for ins_column in ins_columns :

    qry = "SELECT `fullname` FROM `{}`.`{}` WHERE `{}` IS NULL;".format(db_name, tb_name, ins_column)
    cur.execute(qry)
    fullnames = cur.fetchall()

    n = 0

    for fullname in fullnames :
        n += 1
        print('#'*40,
            "\n{2:.01f}%  ({0}/{1}) {3}".format(n, len(fullnames), (n/len(fullnames))*100, os.path.basename(__file__)))
        print ("Starting...   fullname: {}".format(fullname))

        #check file exist...
        if not os.path.exists("{}".format(fullname['fullname'])):
            qry = """DELETE FROM `{0}`.`{1}`  
                    WHERE `{1}`.`fullname` = '{2}';""".format(db_name, tb_name,
                                                        fullname['fullname'])
                
            print("qry: {}".format(qry))
            cur.execute(qry)
            conn.commit()

        else : 
            hdul = fits.open(fullname['fullname'])
            if "{}".format(ins_column) in hdul[0].header :
                qry = """UPDATE `{0}`.`{1}` 
                        SET `{3}`= '{4}'    
                        WHERE `{1}`.`fullname` = '{2}';""".format(db_name, tb_name,
                                                        fullname['fullname'], ins_column, hdul[0].header['{}'.format(ins_column)])
                print("qry: {}".format(qry))
                cur.execute(qry)
                conn.commit()