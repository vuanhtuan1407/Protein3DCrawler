import logging
import urllib.request as req

from Bio import SeqIO
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
# from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.support.ui import WebDriverWait
from tqdm import tqdm

SP6_PATH = './data/train_set.fasta'


def read_fasta(fasta_file=SP6_PATH):
    signal_proteins = []
    records = SeqIO.parse(fasta_file, 'fasta')
    for record in records:
        annotation = str(record.id).split('|')
        seq = str(record.seq)[:len(str(record.seq)) // 2]
        signal_proteins.append({
            "id": annotation[0],
            "organism": annotation[1],
            'label': annotation[2],
            'seq': seq,
        })
    return signal_proteins


if __name__ == '__main__':
    # Read fasta file:
    prots = read_fasta(SP6_PATH)

    # Open driver
    chrome_options = webdriver.ChromeOptions()
    chrome_options.add_argument('--headless=new')
    driver = webdriver.Chrome(options=chrome_options)

    error_ids = []

    # for loop to get file
    for _, prot in tqdm(enumerate(prots)):
        prot_id = prot['id']
        try:
            # access Uniprot page
            driver.get(f'https://www.uniprot.org/uniprotkb/{prot_id}/entry#structure')
            alphafold_url = WebDriverWait(driver, 30).until(
                EC.presence_of_element_located((By.XPATH, "//a[text()='AlphaFold']"))
            ).get_attribute('href')

            # redirect to AlphaFold page
            driver.get(alphafold_url)
            download_url = WebDriverWait(driver, 30).until(
                EC.presence_of_element_located((By.XPATH, "//a[text()='PDB file ']"))
            ).get_attribute('href')

            filename = './data/pdb/' + download_url.split('/')[-1]
            req.urlretrieve(download_url, filename)

        except Exception as e:
            error_ids.append(prot_id)
            logging.exception(e)
