# import time
import json
import os
import time
import urllib.request as req

import numpy as np
import pandas as pd
from Bio import SeqIO
from selenium import webdriver
from selenium.webdriver.common.by import By
# from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.wait import WebDriverWait

from logs.logger import Logger

SP6_PATH = './data/train_set.fasta'
# UNIPROT_URL = 'https://www.uniprot.org/uniprot'
SAVE_DIR = './data/pdb'
# options = webdriver.ChromeOptions()
# driver = webdriver.Chrome(options=options)

logger = Logger(filename='./logs/crawler.log', use_console=True)


def extract_prot_ids_fasta(fasta_file=SP6_PATH):
    prot_ids = []
    records = SeqIO.parse(fasta_file, 'fasta')
    for record in records:
        annotation = str(record.id).split('|')
        prot_ids.append(annotation[0])
    return prot_ids


def load_cache():
    done_ids = np.loadtxt('./done_ids.txt', dtype=str)
    missing_ids = np.loadtxt('./missing_ids.txt', dtype=str)
    return done_ids.tolist(), missing_ids.tolist()


def save_cache(done_ids, missing_ids):
    np.savetxt(fname='done_ids.txt', X=done_ids, fmt='%s')
    np.savetxt(fname='missing_ids.txt', X=missing_ids, fmt='%s')


# First try to download
def download_pdb_by_alphafold():
    prot_ids = extract_prot_ids_fasta()
    if not os.path.exists(SAVE_DIR):
        os.makedirs(SAVE_DIR, exist_ok=True)
    done_ids, missing_ids = load_cache()
    for i, prot_id in enumerate(prot_ids):
        if i > 0 and i % 100 == 0:
            logger.info(f'Processing {i}/{len(prot_ids)} records')
        if prot_id not in done_ids and prot_id not in missing_ids:
            download_url = f'https://alphafold.ebi.ac.uk/files/AF-{prot_id}-F1-model_v4.pdb'
            filename = SAVE_DIR + f'/AF-{prot_id}-F1-model_v4.pdb'
            try:
                req.urlretrieve(download_url, filename)
                logger.info(f'Downloaded at {filename}')
                done_ids.append(prot_id)
            except Exception as e:
                logger.error(f'Failed to download {prot_id}: {e}')
                missing_ids.append(prot_id)

            if i > 0 and i % 100 == 0:
                # logger.info(f'Processing thought {i} records')
                save_cache(done_ids, missing_ids)

        # time.sleep(3)  # prevent bot detection (do not necessary)

    # save after running
    save_cache(done_ids, missing_ids)


# Retry
def retry_download_from_uniprot():
    convert_ids = {}
    done_ids, missing_ids = load_cache()
    options = webdriver.ChromeOptions()
    options.add_argument('--headless')
    driver = webdriver.Chrome(options=options)
    for i, prot_id in enumerate(missing_ids):
        driver.get(f'https://www.uniprot.org/uniprotkb/{prot_id}/entry#structure')
        time.sleep(3)
        try:
            # Retry downloading ids missing from downloading process above by accessing from uniprot page
            current_url = driver.current_url
            af_url = WebDriverWait(driver, 10).until(
                EC.presence_of_element_located((By.XPATH, "//a[text()='AlphaFold']"))
            ).get_attribute('href')
            new_prot_id = current_url.split('/')[-2]  # retrieve new prot_id
            driver.get(af_url)
            time.sleep(3)
            download_url = WebDriverWait(driver, 10).until(
                EC.presence_of_element_located((By.XPATH, "//a[text()='PDB file ']"))
            ).get_attribute('href')
            filename = f'{SAVE_DIR}/AF-{new_prot_id}-F1-model_v4.pdb'
            try:
                req.urlretrieve(download_url, filename)
                missing_ids = missing_ids.remove(prot_id)
                done_ids = done_ids.append(prot_id)
                if prot_id not in convert_ids.keys():
                    convert_ids[prot_id] = new_prot_id
                time.sleep(3)
            except Exception as e:
                logger.error(f'Failed to download {prot_id}: {e}')
        except Exception:
            logger.error('No AlphaFoldl url found')

    # Save convert ids of success retrying
    with open('./convert_ids.json', 'w') as f:
        json.dump(convert_ids, f)


def sync_and_check_existed_files():
    prot_ids = extract_prot_ids_fasta()
    done_ids, missing_ids = load_cache()
    for i, prot_id in enumerate(prot_ids):
        if prot_id not in done_ids:
            if not os.path.exists(f'{SAVE_DIR}/AF-{prot_id}-F1-model_v4.pdb') and prot_id not in missing_ids:
                logger.error(f'Missing protein {prot_id}: Have not processed this protein yet!')
                # missing_ids.append(prot_id)
            elif os.path.exists(f'{SAVE_DIR}/AF-{prot_id}-F1-model_v4.pdb') and prot_id not in missing_ids:
                logger.warning(f'Missing protein {prot_id} in done_ids')
                # done_ids.append(prot_id)
            elif os.path.exists(f'{SAVE_DIR}/AF-{prot_id}-F1-model_v4.pdb') and prot_id in missing_ids:
                logger.warning(f'Missing protein {prot_id} in done_ids but in missing_ids')
                # done_ids.append(prot_id)
                # missing_ids.remove(prot_id)

    # save_cache(done_ids, missing_ids)


def statistic_missing_ids():
    orgs = []
    lbs = []
    partitions = []
    _, missing_ids = load_cache()
    records = SeqIO.parse(SP6_PATH, 'fasta')
    for record in records:
        prot_id, org, lb, partition = str(record.id).split('|')
        if prot_id in missing_ids:
            orgs.append(org)
            lbs.append(lb)
            partitions.append(f'Partition {partition}')

    df = pd.DataFrame({
        "missing_id": missing_ids,
        "org": orgs,
        "lb": lbs,
        "partition": partitions
    })
    df.to_csv('./out/missing_ids.csv', index=True, index_label='index')


if __name__ == '__main__':
    # prot_ids = extract_prot_ids_fasta()
    # download_pdb_by_alphafold()
    # After downloading pdb file by alphafold downloading link, some ids were missing
    # Retry downloading these missing ids from uniprot page
    # retry_download_from_uniprot()
    # sync_and_check_existed_files()
    statistic_missing_ids()
