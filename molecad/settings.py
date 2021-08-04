import logging
import sys
from pathlib import Path

from loguru import logger
from pydantic import BaseSettings, Field, HttpUrl
from pymongo import MongoClient


# from .log import DevelopFormatter, JsonSink


class Settings(BaseSettings):
    env: str = Field("dev", env="ENV")
    project_dir: Path = Field(".", env="PROJ_DIR")
    fetch_dir: Path = Field("./data/fetch", env="FETCH_DIR")
    split_dir: Path = Field("./data/split", env="SPLIT_DIR")

    api_url: HttpUrl = Field("http://127.0.0.1:8000", env="API_URL")
    api_version: str = Field("", env="API_VERSION")
    app_url: HttpUrl = Field("http://127.0.0.1:8050", env="APP_URL")
    app_version: str = Field("", env="APP_VERSION")

    mongo_host: str = Field("127.0.0.1", env="MONGO_HOST")
    mongo_port: int = Field(27017, env="MONGO_PORT")
    mongo_user: str = Field("", env="MONGO_USER")
    mongo_password: str = Field("", env="MONGO_PASSWORD")
    mongo_auth_source: str = Field("admin", env="MONGO_AUTH_SOURCE")

    db_name: str = Field("", env="MONGO_DB_NAME")
    properties: str = Field("properties", env="MONGO_PROPERTIES_COLLECTION")
    molecules: str = Field("molecules", env="MONGO_MOLECULES_COLLECTION")
    mfp_counts: str = Field("mfp_counts", env="MONGO_MFP_COUNTS_COLLECTION")

    component_name: str = "molecad"

    class Config:
        env_file = ".env"
        env_file_encoding = "utf-8"

    @property
    def version(self):
        import importlib_metadata

        return importlib_metadata.version("molecad")

    def setup_logging(self):
        class InterceptHandler(logging.Handler):
            def emit(self, record):
                # Get corresponding Loguru level if it exists
                try:
                    level = logger.level(record.levelname).name
                except ValueError:
                    level = record.levelno

                # Find caller from where originated the logged message
                frame, depth = logging.currentframe(), 2
                while frame.f_code.co_filename == logging.__file__:
                    frame = frame.f_back
                    depth += 1

                logger.opt(depth=depth, exception=record.exc_info).log(level, record.getMessage())

        if self.env == "dev":
            logger.add(
                sys.stdout,
                level="DEBUG",
                format="{name} | {function} | {line}\n {level} | {message}",
                colorize=True,
                serialize=True,
            )
        else:
            logger.add(
                "logs/prod.log",
                level="INFO",
                format="{time} | {name} | {function} | {line} \n{level} | {message} \n{extra[i]}",
                colorize=True,
            )
        # logging.getLogger("uvicorn.access").handlers = [InterceptHandler()]

        logging.basicConfig(handlers=[InterceptHandler()], level=0)

    def get_db(self):
        return MongoClient(
            host=self.mongo_host,
            port=self.mongo_port,
            username=self.mongo_user,
            password=self.mongo_password,
            authSource=self.mongo_auth_source,
        )[self.db_name]

    def get_collections(self):
        db = self.get_db()
        properties = db[self.properties]
        molecules = db[self.molecules]
        mfp_counts = db[self.mfp_counts]
        return properties, molecules, mfp_counts


settings = Settings()
