"""
This file handles management of the riboseqorg database.

The database is a sqlite3 database. The database is created by the
 script. The database is stored in riboseqorg.db  

The database is managed by the class DatabaseManager. This class
    contains methods for creating the database, adding data to the
    database, and querying the database.


"""

import sqlite3
import csv


class DatabaseManager:
    """
    This class handles the database management for the riboseqorg
    """

    def __init__(self, database_name="riboseqorg.db"):
        self.database_name = database_name
        self.conn = sqlite3.connect(self.database_name)
        self.c = self.conn.cursor()
        self.c.execute('SELECT * FROM studies')

# Check if the result set is empty
        if self.c.fetchone() is None:
            self.create_database()
        self.study_fields = """(
                -- The id field is an auto-incrementing primary key that uniquely identifies each sample in the table.
                id INTEGER PRIMARY KEY AUTOINCREMENT, 
                -- The study_id field refers the the corresponding study entry in the studies table
                study_id INTEGER NOT NULL,
                -- The gsm_accession field contains a unique accession number for each sample in the Gene Expression Omnibus (GEO) database.
                gsm_accession TEXT NOT NULL,
                -- The title field contains a brief description of the sample.
                title TEXT NOT NULL,
                -- The organism field specifies the organism from which the sample was derived.
                organism TEXT NOT NULL,
                -- The source field specifies the source of the sample, such as a tissue or cell type.
                source TEXT NOT NULL,
                -- The strain_genotype field specifies the strain or genotype of the organism from which the sample was derived.
                strain_genotype TEXT NOT NULL,
                -- The cell_tissue field specifies the cell or tissue type from which the sample was derived.
                cell_tissue TEXT NOT NULL,
                -- The sample description field contains a detailed description of the sample.
                sample_description TEXT NOT NULL,
                -- The library_strategy field specifies the sequencing library preparation strategy used for the sample.
                library_strategy TEXT NOT NULL,
                -- The ribosome_type field specifies the type of ribosome (e.g. initiating or elongating) targeted in the sample.
                ribosome_type TEXT NOT NULL,
                -- The extraction_protocol field specifies the protocol used to extract RNA from the sample.
                extraction_protocol TEXT NOT NULL,
                -- The tags field can be used to store additional metadata about the sample, such as keywords or tags from GEO
                tags TEXT NOT NULL,
                -- The library_strategy_evidence and ribosome_type_evidence fields contain evidence supporting the values specified in the library_strategy and ribosome_type fields, respectively.
                library_strategy_evidence TEXT NOT NULL,
                ribosome_type_evidence TEXT NOT NULL,
                FOREIGN KEY (study_id) REFERENCES studies(id)
                );
                """
        self.conn.commit()
        self.conn.close()

    def create_database(self):
        """
        Create the database if it does not exist.
        """
        self.c.execute(
            """CREATE TABLE IF NOT EXISTS studies (
                -- A unique identifier for each study
                id INTEGER PRIMARY KEY AUTOINCREMENT, 
                -- GSE Accession number
                gse_accession TEXT NOT NULL, 
                -- Study title
                title TEXT NOT NULL, 
                -- Organism(s) studied (comma-separated)
                organism TEXT NOT NULL, 
                -- Number of samples in the study (as listed in the GSE)
                num_samples INTEGER, 
                -- SRA Accession number
                sra_accession TEXT, 
                -- Date the study was released
                release_date DATE, 
                -- Protocols used in the study
                protocols TEXT, 
                -- Sequencing types used in the study (comma-separated)
                sequencing_types TEXT, 
                -- GSE URL
                gse_url TEXT, 
                -- GSE Supplementary information URL
                gse_supplementary TEXT, 
                -- BioProject number
                bioproject TEXT
                )"""
        )

    def add_study(
        self,
        gse_accession,
        title,
        organism,
        num_samples,
        sra_accession,
        release_date,
        protocols,
        sequencing_types,
        gse,
        gse_supplementary,
        bioproject,
        samples_csv
    ) -> None:
        """
        Add a study to the database
        """
        self.conn = sqlite3.connect(self.database_name)
        self.c = self.conn.cursor()
        #check if study already exists based on gse_accession
        self.c.execute("""SELECT * FROM studies WHERE gse_accession=?""", (gse_accession,))
        if self.c.fetchone() is None:
            print(num_samples)
            self.c.execute(
                """INSERT INTO studies (gse_accession, title, organism, num_samples, sra_accession, release_date, protocols, sequencing_types, gse_url, gse_supplementary, bioproject) VALUES (?,?,?,?,?,?,?,?,?,?,?)""",
                (
                    gse_accession,
                    title,
                    organism,
                    num_samples,
                    sra_accession,
                    release_date,
                    protocols,
                    sequencing_types,
                    gse,
                    gse_supplementary,
                    bioproject,
                ),
            )
            study_id = self.c.lastrowid #get the id of the study that was just added for reference in study table

            self.c.execute(
                """CREATE TABLE IF NOT EXISTS """ + gse_accession + self.study_fields
            )
            with open(samples_csv, "r") as samples_csv:
                csvreader = csv.reader(samples_csv)
                for row in csvreader:
                    self.c.execute(
                        f"""INSERT INTO {gse_accession} (study_id, gsm_accession, title, organism, source, strain_genotype, cell_tissue, sample_description, library_strategy, ribosome_type, extraction_protocol, tags, library_strategy_evidence, ribosome_type_evidence) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
                        (
                            study_id,
                            row[0],
                            row[1],
                            row[2],
                            row[3],
                            row[4],
                            row[5],
                            row[6],
                            row[7],
                            row[8],
                            row[9],
                            row[10],
                            row[11],
                            row[12],
                        ),
                    )
        self.conn.commit()
        self.conn.close()


    def query(self, query, values):
        """
        Query the database
        """
        self.conn = sqlite3.connect(self.database_name)
        self.c = self.conn.cursor()
        self.c.execute(query, values)
        self.conn.commit()
        self.conn.close()
