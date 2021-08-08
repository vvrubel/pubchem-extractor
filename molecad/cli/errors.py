class BadDomainError(Exception):
    pass


class BadNamespaceError(Exception):
    pass


class BadOperationError(Exception):
    pass


class DirExistsError(FileExistsError):
    pass


class EmptySmilesError(Exception):
    pass
