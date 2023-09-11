import ntplib
from datetime import datetime, timezone

def get_ntp_time():

    ntp_pool = ['pool.ntp.org', 'time.nist.gov']

    def call_ntp(serverAddress):
        call = ntplib.NTPClient()
        return call.request(server, version=3)

    for server in ntp_pool:
        response = call_ntp(server)
        print(f"server: {server}")
        print(f"request packet sent (as LOCAL client time, orig_time): {datetime.fromtimestamp(response.orig_time, timezone.utc)}")
        print(f"request packet received (as NTP server time, recv_time): {datetime.fromtimestamp(response.recv_time, timezone.utc)}")
        print(f"response packet sent (as NTP server time, tx_time): {datetime.fromtimestamp(response.tx_time, timezone.utc)}")
        print(f"response packet received (as LOCAL client time, dest_time): {datetime.fromtimestamp(response.dest_time, timezone.utc)}")
        print(f'round trip duration: {response.delay} s')
        print(f'* adjusted time, tx_time + delay/2: {datetime.fromtimestamp(response.tx_time + response.delay/2, timezone.utc)}')
        print(f'* adjusted time, dest_time + offset: {datetime.fromtimestamp(response.dest_time + response.offset, timezone.utc)}')
        print(f'correction to client: {response.delay/2 - response.offset} s\n')
        # for attr in dir(response):
            # if not attr .startswith('_'):
                # print("response.%s = %r" % (attr, getattr(response, attr)))
        print('-')
get_ntp_time()