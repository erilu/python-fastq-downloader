#!/usr/bin/env python3
# source: https://github.com/erilu/python-fastq-downloader/blob/master/fastq_download.py

import os
from glob import glob
import subprocess

# initialization
work_dir = './'
samples = {
    'Het_2': 'SRR2121686',
    'Imm_1': 'SRR2121687',
    'Imm_2': 'SRR2121688',
}

# downloading each given file
for sample_id in samples:
    print('Currently downloading: ' + samples[sample_id])

    # downloading/converting the files
    cmd_prefetch = 'prefetch --output-directory {:s} --progress {:s}'.format(work_dir, samples[sample_id])
    print('\trunning: ' + cmd_prefetch)
    subprocess.call(cmd_prefetch, shell=True)

    cmd_fastqdump = 'fastq-dump --outdir {:s} --skip-technical --readids '.format(work_dir) + \
                    '--read-filter pass --dumpbase --split-3 --clip ' + \
                    '{:s}/{:s}/{:s}.sra'.format(work_dir, samples[sample_id], samples[sample_id])
    print('\trunning: ' + cmd_fastqdump)
    subprocess.call(cmd_fastqdump, shell=True)

    # compressing the fastqs
    for fq_name in glob('{:s}/{:s}*.fastq'.format(work_dir, samples[sample_id])):
        cmd_compress = 'gzip -c {:s} > {:s}/{:s}_{:s}.gz'.format(fq_name, work_dir, sample_id, os.path.basename(fq_name))
        print('\trunning: ' + cmd_compress)
        subprocess.call(cmd_compress, shell=True)
        os.remove(fq_name)

    # clean up
    cmd_rmdir = 'rm -r {:s}/{:s}'.format(work_dir, samples[sample_id])
    print('\trunning: ' + cmd_rmdir)
    subprocess.call(cmd_rmdir, shell=True)

