import pymongo
from pydantic import (
    BaseSettings,
    Field,
)


class Settings(BaseSettings):
    mongo_host = Field("127.0.0.1", env="MONGO_HOST")
    mongo_port = Field(27017, env="MONGO_PORT")
    mongo_user = Field("", env="MONGO_USER")
    mongo_password = Field("", env="MONGO_PASSWORD")
    mongo_auth_source = Field("admin", env="MONGO_AUTH_SOURCE")
    mongo_db_name = Field("molecad", env="MONGO_DB_NAME")
    mongo_collection = Field("", env="MONGO_COLLECTION")
    rdkit_schema = Field("", env="MONGO_SCHEMA_COLLECTION")
    test_collection = Field("test", env="MONGO_TEST_COLLECTION")
    test_schema = Field("test_schema", env="MONGO_TEST_SCHEMA_COLLECTION")

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
        db = pymongo.MongoClient(
            host=self.mongo_host,
            port=self.mongo_port,
            username=self.mongo_user,
            password=self.mongo_password,
            authSource=self.mongo_auth_source
        )[self.mongo_db_name]
        return db

    @property
    def get_prod_collection(self):
        collection = self.get_db[self.mongo_collection]
        schema = self.get_db[self.rdkit_schema]
        return collection, schema

    @property
    def get_dev_collection(self):
        collection = self.get_db["molecules"]
        return collection

    @property
    def get_test_collection(self):
        collection = self.get_db[self.test_collection]
        schema = self.get_db[self.test_schema]
        return collection, schema


setup = Settings()
