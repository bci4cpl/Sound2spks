import os
import csv
import time

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

import multiprocessing

CPU_COUNT = multiprocessing.cpu_count()

np.random.seed(90210)

N_max = 20371
N = 15000
# path = "/home/shay/Sound2spks/WAV_Data/"
path = "/home/shay/TEMP_OFFLINE/WAV_Data/"
seg_file = 'eval_segments_no_header.csv'


def fetch_wav_mock(yt_id, start_time):
    time.sleep(3)
    res = np.random.randint(2)
    if res:
        return yt_id, start_time, f'mock.{yt_id}.wav'
    return yt_id, start_time, None


def fetch_wav(c, yt_id, start_time):
    # https://ytdl-org.github.io/youtube-dl/download.html
    link = os.popen('youtube-dl -i -w --extract-audio --audio-format wav --audio-quality 0 '
                    '--get-url \"https://www.youtube.com/watch?v=' + yt_id + '\"').read()
    if link:
        link = link[:-1]

        dest = path + str(c) + '_' + yt_id + '.wav'
        ret = os.system('ffmpeg -n -ss ' + start_time + ' -i \"' + link + '\" -t 10 ' +
                        ' \"' + dest + '\"')
        if ret in [0, 256]:  # 0: ok, 256: already exists
            return c, yt_id, start_time, dest
    return yt_id, start_time, None


if __name__ == '__main__':
    L = np.random.permutation(N_max)
    L = L[0:N]
    c = 0

    ds = []
    with open(seg_file) as csvfile:
        csvr = csv.reader(csvfile)
        for row in csvr:
            if c in L:
                d = {'c': c, 'yt_id': row[0], 'start_time': row[1]}
                ds.append(d)
            c = c + 1

    df = pd.DataFrame(ds)
    # print(df.head())

    # for i, r in df.iterrows():
    # yt_id, start_time = r['yt_id'], r['start_time']
    # print(yt_id, start_time)
    # print(fetch_wav(yt_id, start_time))

    n_jobs = np.min((N, 5))
    results = Parallel(n_jobs=np.min((n_jobs, CPU_COUNT)), prefer="threads", verbose=10)(
        delayed(fetch_wav)(r['c'], r['yt_id'], r['start_time']) for i, r in df.iterrows())

    print(results)
