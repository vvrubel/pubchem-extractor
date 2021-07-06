from loguru import logger


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


class BadDomainError(BadRequestError):
    def to_dict(self) -> dict:
        return {"error": "Bad domain", "message": str(self)}

    def __str__(self):
        return 'В данной версии сервиса поиск возможен только по базе данных "Compound"'


class BadNamespaceError(BadRequestError):
    def to_dict(self) -> dict:
        return {"error": "Bad namespace", "message": str(self)}

    def __str__(self):
        return 'Пространство имен для поиска по базе данных "Compound" задано некорректно'


class BadOperationError(BadRequestError):
    def to_dict(self) -> dict:
        return {"error": "Bad operation", "message": str(self)}

    def __str__(self):
        return "Ошибка при составлении операции"
