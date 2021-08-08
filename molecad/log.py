import json
import logging
import sys
import traceback as tb
from datetime import datetime

from loguru import logger
from starlette_context import context

from .settings import settings


class LoguruContainer:
    def __init__(self, component_name: str):
        self.component_name = component_name

    @staticmethod
    def add_correlation_id(record):
        try:
            record["extra"].update(
                correlation_id=context.data.get("X-Correlation-ID"),
                request_id=context.data.get("X-Request-ID"),
            )
        except RuntimeError:
            pass


class JsonSink(LoguruContainer):
    def __call__(self, message):
        record = message.record
        self.add_correlation_id(record)
        data = {
            "datetime": record["time"].astimezone().strftime("%Y-%m-%dT%H:%M:%S%z"),
            "timestamp": int(datetime.timestamp(record["time"])),
            "level": record["level"].name,
            "app": self.component_name,
            "msg": record["message"],
            "module": record["module"],
            "function": record["function"],
            "line": record["line"],
            "filename": record["file"].name,
            "logger": record["name"],
        }
        if record["extra"]:
            data["extra"] = record["extra"]
        if record["exception"] is not None:
            data["exception"] = "".join(tb.format_exception(*record["exception"]))

        print(json.dumps(data))


class DevelopFormatter(LoguruContainer):
    def __call__(self, record):
        self.add_correlation_id(record)
        extra = " ".join(f"<lvl>{k}={str(v)}</>" for k, v in record["extra"].items())

        return (
            "<white>{time:YYYY-MM-DD HH:mm:ss.SSS}</> "
            "| <lvl>{level: <8}</> | "
            f"<magenta>{self.component_name}</> "
            "<cyan>{name}</>:<green>{function}</>:<blue>{line}</> "
            "- <lvl>{message}</> - "
            + extra
            + ("\n{exception}\n" if record.get("exception") is not None else "\n")
        )


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


def setup_logging():
    logging.getLogger().handlers = [InterceptHandler()]
    logger.remove()
    if settings.env == "dev":
        develop_fmt = DevelopFormatter(settings.component_name)
        logger.add(sys.stdout, format=develop_fmt, level="DEBUG", backtrace=True, diagnose=True)
    else:
        json_sink = JsonSink(settings.component_name)
        logger.add(json_sink)

    logging.getLogger("uvicorn.access").handlers = [InterceptHandler()]
