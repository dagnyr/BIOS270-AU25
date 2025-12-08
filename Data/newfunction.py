import sqlite3
import pandas as pd
import argparse
import tqdm
import time
import numpy as np
import h5py

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--database_path", type=str, required=True)
    parser.add_argument("--h5_path", type=str, required=True)
    parser.add_argument("--record_id", type=str, required=True)
    parser.add_argument("--output_path", type=str, default="embeddings.npy")
    args = parser.parse_args()
    return args


class BacteriaDatabase:
    def __init__(self, db_path):
        self.db_path = db_path
        self.conn = sqlite3.connect(db_path)

    def get_all_record_ids(self):
        query = "SELECT DISTINCT record_id FROM gff"
        df = self.query(query)
        return df["record_id"].dropna().tolist()

    def get_protein_ids_from_record_id(self, record_id):
        query = "SELECT protein_id FROM gff WHERE record_id = ?"
        df = pd.read_sql(query, self.conn, params=(record_id,))
        return df["protein_id"].dropna().tolist()

    def index_record_ids(self):
        query = "CREATE INDEX IF NOT EXISTS record_id_index ON gff(record_id)"
        self.conn.execute(query)

    def query(self, query):
        return pd.read_sql(query, self.conn)

    def close(self):
        self.conn.close()

    def __del__(self):
        self.close()

def write_protein_ids(protein_ids, output_path):
    with open(output_path, "w") as f:
        f.write("\n".join(protein_ids))
        
def load_embedding_index(hdf5):
    protein_ids = (hdf5["protein_ids"][:]).astype(str)
    return {pid: i for i, pid in enumerate(protein_ids)}

def extract_embeddings(protein_ids, index_map, dataset):
    rows = [index_map[p] for p in protein_ids if p in index_map]
    return dataset[rows, :]

if __name__ == "__main__":
    
    args = parse_args()
    db = BacteriaDatabase(args.database_path)
    
    protein_ids = db.get_protein_ids_from_record_id(args.record_id)
    
    h5f = h5py.File(args.h5_path, "r")
    dataset = h5f["mean_embeddings"]          # only mean
    index_map = load_embedding_index(h5f)
    embeddings = extract_embeddings(protein_ids, index_map, dataset)
    np.save(args.output_path, embeddings)
    h5f.close()

    db.close()
