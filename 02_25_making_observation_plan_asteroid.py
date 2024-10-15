import os
from astropy.io import fits
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from datetime import datetime, timedelta, timezone
import _Python_utilities

#########################################
log_dir = "logs/"
log_file = "{}{}.log".format(log_dir, os.path.basename(__file__)[:-3])
err_log_file = "{}{}_err.log".format(log_dir, os.path.basename(__file__)[:-3])
print ("log_file: {}".format(log_file))
print ("err_log_file: {}".format(err_log_file))
if not os.path.exists('{0}'.format(log_dir)):
    os.makedirs('{0}'.format(log_dir))
#########################################
BASEPATH = Path("./")
save_dir_name = "OBSplan_NEA"
print(f"BASEPATH: {BASEPATH}")


from datetime import datetime, timedelta
start_dt = datetime.strptime("2023-11-13 10:00:00", '%Y-%m-%d %H:%M:%S')
end_dt = datetime.strptime("2023-11-14 10:00:00", '%Y-%m-%d %H:%M:%S')
dt = 1

def datetime_range(start_dt, end_dt, delta):
    '''
        Parameters
        ----------
        start_dt : datetime.datetime
        end_dt : datetime.datetime
        delta : int
    '''
    current = start_dt
    while current < end_dt:
        yield current
        current += delta

dts = [dt for dt in
            datetime_range(start_dt, end_dt,
            timedelta(days=dt))]

print("len(dts):", len(dts))
print("dts", dts)


#############################################
# variables
mpc_code='P64' # Observer Location: GSHS Observatory, Suwon [code: P64]

elev_min=60 #
time_min=15
vmag_min=6
vmag_max=13
list_num=50
sort='trans'

print(len(dts))
print(dts)

df_all = pd.DataFrame()

for dt in dts :
    obs_datetime = dt.strftime('%Y-%m-%dT%H:%M:%S')
    print("obs_datetime :", obs_datetime)
    obs_date = dt.strftime('%Y%m%d')
    print("obs_date :", obs_date)

    if not (BASEPATH/save_dir_name).exists():
        os.mkdir(str(BASEPATH/save_dir_name))
        print (f"{str(BASEPATH/save_dir_name)} is created...")
    else :
        print (f"{str(BASEPATH/save_dir_name)} is already exist...")

    if not (BASEPATH/save_dir_name/obs_date).exists():
        os.mkdir(str(BASEPATH/save_dir_name/obs_date))
        print (f"{str(BASEPATH/save_dir_name/obs_date)} is created...")
    else :
        print (f"{str(BASEPATH/save_dir_name/obs_date)} is already exist...")
        

    # Define API URL and SPK filename:
    url = 'https://ssd-api.jpl.nasa.gov/sbwobs.api'
    spk_filename = 'spk_file.bsp'

    # Get the requested SPK-ID from the command-line:
    if (len(sys.argv)) == 1:
        print("please specify SPK-ID on the command-line")
        sys.exit(2)
    spkid = sys.argv[1]

    # Build the appropriate URL for this API request:
    # IMPORTANT: You must encode the "=" as "%3D" and the ";" as "%3B" in the
    #            Horizons COMMAND parameter specification.
    url += f"?sb-kind=a&mpc-code={mpc_code}&obs-time={str(obs_datetime)}&elev-min={elev_min}&time-min={time_min}"
    url += f"&vmag-max={vmag_max}&vmag-min={vmag_min}&optical=true&fmt-ra-dec=true&mag-required=true&output-sort={sort}&maxoutput={list_num}"

    print("url :", url)

    import json
    import requests
    # Submit the API request and decode the JSON-response:
    response = requests.get(url)
    try:
        data = json.loads(response.text)
        print(data)
    except ValueError:
        print("Unable to decode JSON results")

    df = pd.DataFrame.from_dict(data['data'])
    df = df.set_axis((data['fields']), axis=1)
    df['OBSdate(UT)'] = obs_date
    #print(df)

    df.to_csv(str(BASEPATH/save_dir_name/f"OBSPlan_NEA_{dt.strftime('%Y%m%d')}.csv"))

    df_all = pd.concat([df_all, df], axis = 0)

df_all.reset_index(inplace=True)
print("df_all :", df_all)

from pathlib import Path
from astroquery.mpc import MPC
from pprint import pprint
import re

for i, row  in df.iterrows():
    print("row['Full name'] :", row["Full name"])

    #asteroid_names = re.findall('\(([^)]+)', df['Full name'][0])   #extracts string in bracket()
    #print (asteroid_names[0])
    asteroid_id = df['Full name'][i].split(" ")[0]   #extracts string
    print(asteroid_id)

    eph = MPC.get_ephemeris(asteroid_id, step='1h', number=30 * 24)
    print(f"Proper motion: {eph['Proper motion'].quantity.to('arcmin/h').max():.03f}")
    eph = MPC.get_ephemeris(asteroid_id, location='P64',
                        ra_format={'sep': ':', 'unit': 'hourangle', 'precision': 1},
                        dec_format={'sep': ':', 'precision': 0},
                        step='1h')
    print("type(eph): ", type(eph))
    #print("eph: ", eph)

    df_eph = eph.to_pandas()
    df_eph['Asteroid ID'] = asteroid_id
    df_eph['Asteroid name'] = df['Full name'][i]
    df_eph.to_csv(str(BASEPATH/save_dir_name/ f"eph_asteroid_{asteroid_id}.csv"))
