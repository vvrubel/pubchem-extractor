from pydantic import BaseSettings, Field
from pymongo import MongoClient


class Settings(BaseSettings):
    mongo_host = Field("127.0.0.1", env="MONGO_HOST")
    mongo_port = Field(27017, env="MONGO_PORT")
    mongo_user: str = Field("", env="MONGO_USER")
    mongo_password: str = Field("", env="MONGO_PASSWORD")
    mongo_auth_source = Field("admin", env="MONGO_AUTH_SOURCE")

    db_name: str = Field("", env="MONGO_DB_NAME")
    pubchem = Field("pubchem", env="MONGO_PUBCHEM_COLLECTION")
    molecules = Field("molecules", env="MONGO_MOLECULES_COLLECTION")
    mfp_counts = Field("mfp_counts", env="MONGO_MFP_COUNTS_COLLECTION")
    permutations = Field("permutations", env="MONGO_PERMUTATIONS_COLLECTION")

    project_dir = Field(".", env="PROJ_DIR")
    fetch_dir = Field("./data/fetch", env="FETCH_DIR")
    split_dir = Field("./data/split", env="SPLIT_DIR")

    component_name = "molecad"
    api_version: str = "/v0"

    class Config:
        env_file = ".env"
        env_file_encoding = "utf-8"

    @property
    def get_db(self):
        return MongoClient(
            host=self.mongo_host,
            port=self.mongo_port,
            username=self.mongo_user,
            password=self.mongo_password,
            authSource=self.mongo_auth_source,
        )[self.db_name]

    @property
    def get_collections(self):
        db = self.get_db
        pubchem = db[self.pubchem]
        molecules = db[self.molecules]
        mfp_counts = db[self.mfp_counts]
        permutations = db[self.permutations]
        return pubchem, molecules, mfp_counts, permutations


settings = Settings()
