import json
import sys
import traceback as tb
from datetime import datetime

from loguru import logger
from starlette_context import context


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
            "<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</> "
            "| <lvl>{level: <8}</> | "
            f"<cyan>{self.component_name}</> "
            "<cyan>{name}</>:<cyan>{function}</>:<cyan>{line}</> "
            "- <lvl>{message}</> - "
            + extra
            + ("\n{exception}\n" if record.get("exception") is not None else "\n")
        )


#
# logger.remove()
# logger.add(
#     sys.stdout,
#     colorize=True,
#     format="<green>{time:HH:mm:ss}</green> | {level} | <level>{message}</level>",
# )
# logger.add("logs/fetch.log", level="DEBUG", format="{time} {level} {message}")

logger.configure(
    handlers=[
        dict(sink=sys.stderr, format="[{time}] | {level} | {message}"),
        dict(sink="file.log", enqueue=True, serialize=True),
    ],
    levels=[dict(name="NEW", no=13, icon="¤", color="")],
    extra={"common_to_all": "default"},
    patcher=lambda record: record["extra"].update(some_value=42),
    activation=[("my_module.secret", False), ("another_library.module", True)],
)

# Set a default "extra" dict to logger across all modules, without "bind()"
extra = {"context": "foo"}
logger.configure(extra=extra)
logger.add(sys.stderr, format="{extra[context]} - {message}")
logger.info("Context without bind")
# => "foo - Context without bind"
logger.bind(context="bar").info("Suppress global context")
# => "bar - Suppress global context"

fmt = "{time} - {name} - {level} - {message}"
logger.add("spam.log", level="DEBUG", format=fmt)
logger.add(sys.stderr, level="ERROR", format=fmt)
