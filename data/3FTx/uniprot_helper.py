import json
import sys
import time
from io import StringIO
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

# Documentation: https://www.uniprot.org/help/api
WEBSITE_API = "https://rest.uniprot.org"

# Documentation: https://www.ebi.ac.uk/proteins/api/doc/
PROTEINS_API = "https://www.ebi.ac.uk/proteins/api"

# Documentation here https://www.ebi.ac.uk/seqdb/confluence/pages/viewpage.action?pageId=94147939#NCBIBLAST+HelpandDocumentation-RESTAPI
BLAST_API = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast"

TAXA_API = "https://www.ebi.ac.uk/proteins/api/taxonomy"

# Helper function to download data
def get_url(url, **kwargs):
    response = requests.get(url, **kwargs)

    if not response.ok:
        print(response.text)
        response.raise_for_status()
        sys.exit()

    return response


def get_uniprot_entry(acc_id: str, format_type: str = "json"):
    """Get the result of a single UniProt entry by accession ID"""
    url = f"{WEBSITE_API}/uniprotkb/{acc_id}.{format_type}"
    response = get_url(url=url)
    return response


def get_tax_id(taxa: str) -> str:
    """For a given taxa name, search and return its taxa ID.

    Return the taxa ID as an integer if it exists, else return None.
    """
    url = f"{TAXA_API}/name/{taxa}"
    response = get_url(url)
    try:
        taxonomies = response.json()["taxonomies"]
    except requests.HTTPError as err:
        txt = err
        if str(txt).startswith("404 Client Error: Not Found for url"):
            taxonomies = None
        else:
            raise Exception(err)

    if taxonomies is None:
        tax_id = None
    elif len(taxonomies) == 1:
        tax_id = taxonomies[0]["taxonomyId"]
    else:
        raise Exception(f"{taxa} found more than one taxa.")
    return tax_id


def post_blast_job(seq, taxa_id):
    """Posts a BLAST request"""
    respond = requests.post(
        f"{BLAST_API}/run",
        data=dict(
            alignments=5,
            database="uniprotkb",
            email="tobias.senoner@tum.de",
            # exp="1e-10",
            filter="F",
            gapalign="false",
            matrix="BLOSUM62",
            program="blastp",
            scores=5,
            sequence=seq,
            stype="protein",
            taxids=taxa_id,
        ),
    )
    job_id = respond.text
    return job_id


def get_job_status(job_id):
    """Get the job status. (RUNNING, FINISHED)"""
    response = get_url(f"{BLAST_API}/status/{job_id}")
    return response.text


def get_job_result(job_id: str, result_type: str = "out"):
    """Get the results returned by a finished job by job_id

    Args:
        job_id (str): job identifier returned by the POST request
        resulttypes (str): one of -> out, json, ids, accs, xml, taxids, tsv,
                           error, sequence, visual-svg,complete-visual-svg,
                           visual-png, complete-visual-png

    Returns:
        response (str): response of the requested `job_id`
    """
    return get_url(f"{BLAST_API}/result/{job_id}/{result_type}")


class BatchBlaster:
    """Submit <= 30 jobs at a time. Wait for results before next batch."""

    def __init__(self, entries: list[dict], out_dir: Path) -> None:
        # entries = list of dict with keys: seq, taxon_id, fasta_id
        # TODO: assert that entries are in the correct format (have corect/needed keys)
        self.entries = entries
        self.out_dir = out_dir
        self.entries_done = [file.stem for file in self.out_dir.glob("*.json")]

        self.batch_size = 30
        self.timeout_time = 300
        self.batch_jobs = list()

    @staticmethod
    def batch(iterable, size=1):
        iter_len = len(iterable)
        for idx in range(0, iter_len, size):
            yield iterable[idx : min(idx + size, iter_len)]

    def _save_result(self, job_id: str, entry_id: str):
        job_result = get_job_result(job_id=job_id, result_type="json")
        json_path = self.out_dir / f"{entry_id}.json"
        with open(json_path, "w") as json_handler:
            json.dump(job_result.json(), fp=json_handler, indent=4)

    def _wait_for_jobs2finish(self, pbar):
        nr_done = 0
        seconds_no_update = 0
        while nr_done != (len(self.jobs)):
            start_done = nr_done
            time.sleep(60)
            for entry_id, job in self.jobs.items():
                if job["status"] != "FINISHED":
                    status = get_job_status(job_id=job["job_id"])
                    if status == "FINISHED":
                        self.jobs[entry_id]["status"] = status
                        self._save_result(
                            job_id=job["job_id"], entry_id=entry_id
                        )
                        self.entries_done.append(entry_id)
                        nr_done += 1
                        pbar.update(1)
            if start_done == nr_done:
                seconds_no_update += 60
            # if not a single job got updated in the timeout time exit loop
            if seconds_no_update >= self.timeout_time:
                print("Timeout")
                break

    def run_batch(self):
        pbar = tqdm(total=len(self.entries), desc="BLAST")
        for batch_entries in self.batch(self.entries, size=self.batch_size):
            # submit a `batch_size` of jobs at a time
            jobs = {}
            for entry in batch_entries:
                entry_id = entry["fasta_id"]
                if entry_id in self.entries_done:
                    pbar.update(1)
                    continue
                else:
                    job_id = post_blast_job(
                        seq=entry["seq"], taxa_id=entry["taxon_id"]
                    )
                    jobs[entry_id] = dict(job_id=job_id, status="RUNNING")
            self.jobs = jobs

            self._wait_for_jobs2finish(pbar=pbar)
        pbar.close()
