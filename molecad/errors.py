import pymongo.errors
from loguru import logger


class BadDomainError(Exception):
    pass


class BadNamespaceError(Exception):
    pass


class BadOperationError(Exception):
    pass


class DirExistsError(FileExistsError):
    pass


class CreateIndexError(Exception):
    pass


class BaseAppException(Exception):
    def __init__(self, message: str = "") -> None:
        super().__init__(message)
        self.message = message
        logger.trace("Exception '{}' was raised. {}.", self.__class__.__name__, self.message)

    @property
    def error_code(self) -> int:
        raise NotImplementedError("Error code should be implemented.")

    def to_dict(self) -> dict:
        raise NotImplementedError("`to_dict` is not implemented.")


class BadRequestError(BaseAppException):
    def to_dict(self) -> dict:
        return {"error": "Bad request"}

    @property
    def error_code(self) -> int:
        return 400
