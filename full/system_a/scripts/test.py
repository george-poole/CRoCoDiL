import datetime
import time
import random

t = random.uniform(0, 5)
print(f'sleeping {t}')
time.sleep(t)
print(str(datetime.datetime.now()), flush=True)