
import enum

__all__ = [
    "PerturbationTheoryException",
    "Settings",
    "LogLevels",
    "MemoryConstraints"
]
class PerturbationTheoryException(Exception):
    pass


class Settings:

    non_zero_cutoff = 1.0e-14


class LogLevels(enum.Enum):
    Data = 1

class MemoryConstraints(enum.Enum):
    ...