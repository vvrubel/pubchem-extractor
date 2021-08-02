import logging
import sys
from pathlib import Path

from loguru import logger
from pydantic import BaseSettings, Field, HttpUrl

from .log import DevelopFormatter, JsonSink


class Settings(BaseSettings):
    env: str = Field("PROD", env="ENV")
    project_dir: Path = Field(".", env="PROJ_DIR")
    fetch_dir: Path = Field("./data/fetch", env="FETCH_DIR")
    split_dir: Path = Field("./data/split", env="SPLIT_DIR")
    log_cli: Path = Field("./tmp/timer.log", env="LOG_CLI")

    api_url: HttpUrl = Field("http://127.0.0.1:8000", env="API_URL")
    api_version: str = Field("", env="API_VERSION")
    app_url: HttpUrl = Field("http://127.0.0.1:8050", env="APP_URL")

    mongo_host: str = Field("127.0.0.1", env="MONGO_HOST")
    mongo_port: int = Field(27017, env="MONGO_PORT")
    mongo_user: str = Field("", env="MONGO_USER")
    mongo_password: str = Field("", env="MONGO_PASSWORD")
    mongo_auth_source: str = Field("admin", env="MONGO_AUTH_SOURCE")

    db_name: str = Field("", env="MONGO_DB_NAME")
    prop_collection: str = Field("properties", env="MONGO_PROPERTIES_COLLECTION")
    mol_collection: str = Field("molecules", env="MONGO_MOLECULES_COLLECTION")
    mfp_collection: str = Field("mfp_counts", env="MONGO_MFP_COUNTS_COLLECTION")

    component_name: str = "molecad"

    class Config:
        env_file = ".env"
        env_file_encoding = "utf-8"

    @property
    def mongo_uri(self):
        return (
            f"mongodb://{self.mongo_user}:{self.mongo_password}@{self.mongo_host}:"
            f"{self.mongo_port}/{self.mongo_auth_source}"
        )

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

        logging.getLogger().handlers = [InterceptHandler()]
        logger.remove()
        if self.env == "DEV":
            develop_fmt = DevelopFormatter(self.component_name)
            logger.add(sys.stdout, format=develop_fmt, level="DEBUG", backtrace=True, diagnose=True)
        else:
            json_sink = JsonSink(self.component_name)
            logger.add(json_sink)

        logging.getLogger("uvicorn.access").handlers = [InterceptHandler()]
