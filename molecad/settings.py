import pymongo
from pydantic import (
    BaseSettings,
    Field,
)


class Settings(BaseSettings):
    env: str = Field("production", env="ENV")
    api_version: str = "/v0"
    service_url: str = Field("http://localhost:8888", env="SERVICE_URL")

    mongo_host = Field("127.0.0.1", env="MONGO_HOST")
    mongo_port = Field(27017, env="MONGO_PORT")
    mongo_user = Field("", env="MONGO_USER")
    mongo_password = Field("", env="MONGO_PASSWORD")
    mongo_db_name = Field("molecad", env="MONGO_DB_NAME")
    mongo_db_collection = Field("molecules", env="MONGO_DB_COLLECTION")
    mongo_auth_source = Field("admin", env="MONGO_AUTH_SOURCE")

    project_dir = Field("", env="PROJ_DIR")
    data_dir = Field("", env="OUT_DIR")
    saved_json = Field("", env="JSON")

    component_name = "molecad"

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
            authSource=self.mongo_db_name
        )[self.mongo_db_name]
        return db


setup = Settings()
