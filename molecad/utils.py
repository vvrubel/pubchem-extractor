import functools
import time
import urllib.parse
from typing import Callable, Dict, TypeVar, Union

from loguru import logger

T = TypeVar("T")


def timer(func: Callable[..., T]) -> Callable[..., T]:
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        logger.start()
        start = time.monotonic()
        val = func(*args, **kwargs)
        end = time.monotonic()
        work_time = end - start
        hours = work_time // 3600
        minutes = (work_time % 3600) // 60
        seconds = (work_time % 3600) % 60
        logger.success(
            f"Время выполнения {func.__name__!r}:\n "
            f"{hours} часов, {minutes} минут, {seconds:.2f} секунд."
        )
        return val

    return wrapper


def logger_wraps(*, entry=True, exit=True, level="DEBUG"):
    def wrapper(func):
        name = func.__name__

        @functools.wraps(func)
        def wrapped(*args, **kwargs):
            logger_ = logger.opt(depth=1)
            if entry:
                logger_.log(level, "Entering '{}' (args={}, kwargs={})", name, args, kwargs)
            result = func(*args, **kwargs)
            if exit:
                logger_.log(level, "Exiting '{}' (result={})", name, result)
            return result

        return wrapped

    return wrapper


def url_encoder(route: str, query: Dict[str, Union[str, int]]) -> str:
    """
    Энкодер URL-адреса для api. Нужен для того, чтобы не использовать сторонние сервисы для
    генерирования URL.
    :param route: Путь до ручки api, например: '/v1/compound'.
    :param query: Параметры, передаваемые в URL-адрес в формате словаря.
    Например: ``query = {'smiles': 'NC(=O)N', 'skip': 0, 'limit': 10}``.
    :return: Строка URL-адреса.
    """
    params = urllib.parse.urlencode(query)
    return f"http://127.0.0.1:8000{route}?{params}"
