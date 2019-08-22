import os
import csv
import numpy as np

N_max = 20371
N = 10
path = ""
seg_file = 'eval_segments_no_header.csv'

L = np.random.permutation(N_max)
L = L[0:N]
c = 0
with open(seg_file) as csvfile:
    csvr = csv.reader(csvfile)
    for row in csvr:
        if c in L:
            yt_id = row[0]
            start_time = row[1]
            link = os.popen('youtube-dl -i -w --extract-audio --audio-format wav --audio-quality 0 '
                            '--get-url \"https://www.youtube.com/watch?v=' + yt_id + '\"').read()
            link = link[0:-1]
            os.system('ffmpeg -n -ss ' + start_time + ' -i \"' + link + '\" -t 10 ' +
                      ' \"' + path + str(c) + '_' + yt_id + '.wav\"')
        c = c + 1
